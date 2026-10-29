"""
Tests for :mod:`nmrtools`.

Two groups:

* **Unit tests** on synthetic data with analytically known answers. They
  need no dataset and always run.
* **Data tests** that load the real datasets of the Phase 0 regression
  cases and check ``nmrtools`` against the reference results of the legacy
  ``Datos`` classes (``tests/golden``). Where the only intended difference
  is the ppm axis (nmrglue's instead of the legacy one), the test rebuilds
  the legacy axis with :func:`legacy_ppm_axis` so that the analysis
  functions themselves can be compared to floating-point precision.
"""
import dataclasses
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import nmrtools as nt  # noqa: E402
from make_golden import golden_file  # noqa: E402
from regression_cases import DENDRITES, IMEC_OCV, PEO, RC1, data_path  # noqa: E402

RTOL = 1e-9


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def golden(name):
    """
    Load a Phase 0 reference file as a plain dictionary.

    Parameters
    ----------
    name : str
        Case name in ``tests/golden``.

    Returns
    -------
    dict of str to numpy.ndarray
    """
    with np.load(golden_file(name)) as g:
        return {k: g[k] for k in g.files}


def require(*relative):
    """
    Skip the calling test when a dataset folder is not on this machine.

    Parameters
    ----------
    *relative : str
        Dataset folders relative to ``regression_cases.DATA_ROOT``.
    """
    for rel in relative:
        if not os.path.isdir(data_path(rel)):
            pytest.skip(f"dataset not on this machine: {rel}")


def legacy_ppm_axis(dic):
    """
    Rebuild the direct ppm axis exactly as ``DatosProcesados.crear_ppmAxis``.

    Parameters
    ----------
    dic : dict
        nmrglue dictionary of a :class:`nmrtools.Dataset`.

    Returns
    -------
    numpy.ndarray
        ``linspace(OFFSET, OFFSET - SW, FTSIZE)``, with SW (ppm) from acqus.
    """
    offset = dic["procs"]["OFFSET"]
    return np.linspace(offset, offset - dic["acqus"]["SW"], dic["procs"]["FTSIZE"])


def assert_summary_matches(arr, ref, prefix):
    """
    Check an array against a ``regression_cases.summarize`` record.

    Parameters
    ----------
    arr : numpy.ndarray
        Freshly computed array.
    ref : dict
        Golden dictionary that contains the ``prefix.*`` entries.
    prefix : str
        Name the array was summarised under.
    """
    from regression_cases import summarize
    out = {}
    summarize(out, prefix, arr)
    for key, value in out.items():
        expected = ref[key]
        scale = float(np.max(np.abs(expected))) if expected.size else 0.0
        np.testing.assert_allclose(value, expected, rtol=RTOL, atol=RTOL * scale,
                                   err_msg=f"'{key}' differs")


def lorentzian(x, x0, width, area):
    """Lorentzian line with a given area and full width at half maximum."""
    hw = width / 2
    return area / np.pi * hw / ((x - x0) ** 2 + hw ** 2)


# ---------------------------------------------------------------------------
# Unit tests: parameters
# ---------------------------------------------------------------------------
DIC = {
    "acqus": {"NS": 16, "RG": 101.0, "D": [0.0, 2.5, 0.0], "SFO1": 155.5,
              "NUC1": "<7Li>", "PULPROG": "<zg30>", "GPZ": [0, 0, 0, 0, 0, 0, 80.0]},
    "procs": {"SI": 1024, "OFFSET": 100.0},
    "acqu2s": {"TD": 12},
}


def test_param_direct_and_indexed():
    """Plain names are looked up directly, ``D1``-style names index arrays."""
    assert nt._param(DIC, "NS") == 16
    assert nt._param(DIC, "SFO1") == 155.5      # key ending in a digit, not an index
    assert nt._param(DIC, "D1") == 2.5
    assert nt._param(DIC, "GPZ6") == 80.0
    assert nt._param(DIC, "SI") == 1024         # falls through to procs
    assert nt._param(DIC, "TD", dim=2) == 12


def test_param_missing_raises():
    """Unknown parameters raise a KeyError naming the files searched."""
    with pytest.raises(KeyError, match="acqus"):
        nt._param(DIC, "NOPE")


def test_gamma_labels():
    """Bruker labels with angle brackets and any case are accepted."""
    assert nt.gamma("<7Li>") == nt.gamma("7li") == 103.962e6
    with pytest.raises(KeyError):
        nt.gamma("999Xx")


def test_gamma_matches_legacy_table():
    """The table reproduces every value of the legacy ``Datos.set_gamma``."""
    nuclei = ["1H", "6Li", "7Li", "13C", "19F", "23Na", "31P", "29Si", "79Br"]
    np.testing.assert_allclose([nt.gamma(n) for n in nuclei],
                               golden("synthetic_gamma")["gamma"], rtol=1e-15)


# ---------------------------------------------------------------------------
# Unit tests: spectral functions
# ---------------------------------------------------------------------------
PPM = np.linspace(300, -100, 40001)   # descending, like Bruker


def test_crop_inclusive_and_copies():
    """Both limits are kept, the input is untouched, order does not matter."""
    spec = np.arange(PPM.size, dtype=float)
    ppm, cut = nt.crop(PPM, spec, (0, 10))
    assert ppm[0] == 10 and ppm[-1] == 0
    np.testing.assert_array_equal(nt.crop(PPM, spec, (10, 0))[1], cut)
    cut[:] = -1
    assert spec.min() == 0


def test_crop_2d_and_empty_range():
    """Cropping acts on the last axis; an empty range is an error."""
    spec = np.ones((3, PPM.size))
    assert nt.crop(PPM, spec, (0, 10))[1].shape == (3, 1001)
    with pytest.raises(ValueError):
        nt.crop(PPM, spec, (1000, 2000))


@pytest.mark.parametrize("axis", ["descending", "ascending"])
def test_integrate_known_area(axis):
    """A Lorentzian of area 5 integrates to 5, whatever the axis direction."""
    ppm = PPM if axis == "descending" else PPM[::-1]
    spec = lorentzian(ppm, 100.0, 0.5, 5.0)
    area = nt.integrate(ppm, spec, (50, 150))
    assert area == pytest.approx(5.0 * (2 / np.pi) * np.arctan(50 / 0.25), rel=1e-6)


def test_integrate_rows_and_parts():
    """2D input gives one integral per row; parts select the component."""
    line = lorentzian(PPM, 100.0, 0.5, 1.0)
    spec = np.vstack([line, 2 * line]) * (1 + 1j)
    real = nt.integrate(PPM, spec, (50, 150))
    np.testing.assert_allclose(real[1] / real[0], 2.0)
    np.testing.assert_allclose(nt.integrate(PPM, spec, (50, 150), "imag"), real)
    np.testing.assert_allclose(nt.integrate(PPM, spec, (50, 150), "abs"), np.sqrt(2) * real)
    with pytest.raises(ValueError):
        nt.integrate(PPM, spec, part="phase")


def test_mean_ppm_symmetric_peak():
    """The centre of mass of a symmetric peak is its position."""
    spec = lorentzian(PPM, 42.0, 1.0, 1.0)
    assert nt.mean_ppm(PPM, spec, (32, 52)) == pytest.approx(42.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Unit tests: delay lists and diffusion
# ---------------------------------------------------------------------------
def test_read_vdlist_units(tmp_path):
    """Plain numbers are seconds; u/m/s suffixes are converted; blanks skipped."""
    (tmp_path / "vdlist").write_text("0.001\n10m\n500u\n\n2s\n")
    np.testing.assert_allclose(nt.read_vdlist(tmp_path), [1e-3, 1e-2, 5e-4, 2.0])


def test_vd_factor(tmp_path):
    """``vd*N`` in the pulse program is found; otherwise the default is used."""
    (tmp_path / "pulseprogram").write_text("; satrec\n  vd*2\n go=2\n")
    assert nt.vd_factor(tmp_path) == 2.0
    (tmp_path / "pulseprogram").write_text("; zg\n d1\n go=2\n")
    assert nt.vd_factor(tmp_path) == 1.0


@pytest.mark.parametrize("shape", ["sin", "rect"])
@pytest.mark.parametrize("factor_b", [1.0, 1.3])
def test_bvalues_match_legacy(shape, factor_b):
    """b-values equal the legacy ones (which were in 1e9 s/m^2, times factor_b)."""
    ref = golden("synthetic_bvalues")
    gp = np.linspace(0, 80.0, 16)
    b = nt.bvalues(gp, delta=2.0e-3, big_delta=50.0e-3, gamma=103.962e6, shape=shape)
    np.testing.assert_allclose(b * 1e-9 * factor_b, ref[f"{shape}.{factor_b}.bvalue"], rtol=1e-12)


def test_ste_bvalues_reads_parameters():
    """δ = 2·P30, Δ = D20, ramp 0..GPZ6 over TD(F1) points."""
    dic = {"acqus": {"P": [0.0] * 30 + [1000.0], "D": [0.0] * 20 + [0.05],
                     "GPZ": [0, 0, 0, 0, 0, 0, 80.0], "NUC1": "<7Li>"},
           "acqu2s": {"TD": 16}}

    class Holder(nt._ParamMixin):
        """Minimal object with a ``dic`` attribute, as Dataset/Fid have."""

        def __init__(self, d):
            self.dic = d

    gp, b = nt.ste_bvalues(Holder(dic))
    np.testing.assert_allclose(gp, np.linspace(0, 80, 16))
    np.testing.assert_allclose(b, nt.bvalues(gp, 2e-3, 0.05, 103.962e6))


# ---------------------------------------------------------------------------
# Unit tests: fitting
# ---------------------------------------------------------------------------
TAU = np.geomspace(1e-3, 10, 25)


def test_exp_fit_mono_recovery():
    """Saturation recovery with T1 = 0.5 s is recovered exactly from clean data."""
    y = 100 * (1 - np.exp(-TAU / 0.5))
    fit = nt.exp_fit(TAU, y, "mono", "recovery")
    assert fit.params["tau"] == pytest.approx(0.5, rel=1e-6)
    assert fit.params["A"] == pytest.approx(-100, rel=1e-6)
    assert fit.r2 == pytest.approx(1.0)
    assert set(fit.ci95) == {"A", "tau", "y0"}
    assert np.max(np.abs(fit.residuals)) < 1e-6


def test_exp_fit_mono_decay_unsorted():
    """Decays are fitted with A >= 0, and unsorted input is handled."""
    rng = np.random.default_rng(0)
    order = rng.permutation(TAU.size)
    y = 3 * np.exp(-TAU / 0.02) + 0.1
    fit = nt.exp_fit(TAU[order], y[order], "mono", "decay")
    assert fit.params["tau"] == pytest.approx(0.02, rel=1e-6)
    assert np.all(np.diff(fit.x) > 0)


def test_exp_fit_bi_recovery():
    """Two well-separated components are both recovered."""
    y = 1000 - 600 * np.exp(-TAU / 0.01) - 400 * np.exp(-TAU / 1.0)
    fit = nt.exp_fit(TAU, y, "bi", "recovery")
    taus = sorted([fit.params["tau1"], fit.params["tau2"]])
    np.testing.assert_allclose(taus, [0.01, 1.0], rtol=1e-4)


def test_fit_curve_is_log_spaced_for_wide_ranges():
    """``curve()`` spans the data and uses log spacing over many decades."""
    fit = nt.exp_fit(TAU, 1 - np.exp(-TAU / 0.5))
    x, y = fit.curve(n=50)
    assert x[0] == pytest.approx(TAU[0]) and x[-1] == pytest.approx(TAU[-1])
    np.testing.assert_allclose(np.diff(np.log(x)), np.log(x[1] / x[0]))
    np.testing.assert_allclose(y, fit.predict(x))


def test_diffusion_fit():
    """D and S0 of a clean Stejskal–Tanner decay are recovered."""
    b = np.linspace(0, 5e9, 12)
    fit = nt.diffusion_fit(b, 7.0 * np.exp(-b * 1.2e-10))
    assert fit.params["D"] == pytest.approx(1.2e-10, rel=1e-9)
    assert fit.params["S0"] == pytest.approx(7.0, rel=1e-9)


# ---------------------------------------------------------------------------
# Data tests
# ---------------------------------------------------------------------------
@pytest.mark.data
@pytest.mark.parametrize("name, rel, procno", [
    ("rc1_pseudo2d", RC1, 1),
    ("rc1_1d_procno22", RC1, 22),
    ("peo_t1", f"{PEO}/20", 1),
    ("peo_1d_19f", f"{PEO}/27", 1),
    ("imec_ocv_pseudo2d", IMEC_OCV, 1),
])
def test_load_spectrum_matches_legacy(name, rel, procno):
    """Spectra and NS·RG normalisation are identical to the legacy reader."""
    require(rel)
    ds = nt.load(data_path(rel), procno=procno)
    assert_summary_matches(ds.spec, golden(name), "spec")


@pytest.mark.data
@pytest.mark.parametrize("rel, procno", [(RC1, 1), (RC1, 22), (f"{PEO}/27", 1)])
def test_ppm_axis_is_bruker_definition(rel, procno):
    """The ppm axis is ``OFFSET - i * SW_p / (SF * SI)``."""
    require(rel)
    ds = nt.load(data_path(rel), procno=procno)
    pr = ds.dic["procs"]
    i = np.arange(ds.ppm.size)
    np.testing.assert_allclose(ds.ppm, pr["OFFSET"] - i * pr["SW_p"] / pr["SF"] / pr["SI"],
                               rtol=0, atol=1e-9)


@pytest.mark.data
def test_dataset_parameters_and_metadata():
    """Parameter shortcuts and metadata of a real dataset."""
    require(RC1)
    ds = nt.load(data_path(RC1), procno=1)
    assert ds.ndim == 2 and ds.ppm_ind is not None
    assert ds.p("NS") == golden("rc1_pseudo2d")["acqus.NS"]
    assert ds.p("D1") == golden("rc1_pseudo2d")["acqus.D1"]
    assert ds.nucleus == "7Li" and ds.pulprog
    assert ds.norm == ds.p("NS") * ds.p("RG")


@pytest.mark.data
def test_dataset_is_immutable():
    """Cropping returns a new dataset; fields cannot be reassigned."""
    require(RC1)
    ds = nt.load(data_path(RC1), procno=22)
    small = ds.crop((300, 200))
    assert ds.spec.shape != small.spec.shape
    assert small.spec.shape[-1] == small.ppm.size
    with pytest.raises(dataclasses.FrozenInstanceError):
        ds.spec = small.spec


@pytest.mark.data
def test_integration_matches_legacy_on_legacy_axis():
    """With the legacy axis, integrals reproduce the legacy numbers exactly."""
    require(RC1)
    ref2d, ref1d = golden("rc1_pseudo2d"), golden("rc1_1d_procno22")
    ds2 = nt.load(data_path(RC1), procno=1)
    ax = legacy_ppm_axis(ds2.dic)
    for part, key in (("real", "int_metal"), ("abs", "int_metal_abs"), ("imag", "int_metal_imag")):
        np.testing.assert_allclose(nt.integrate(ax, ds2.spec, (280, 220), part), ref2d[key], rtol=RTOL)
    ds1 = nt.load(data_path(RC1), procno=22)
    ax = legacy_ppm_axis(ds1.dic)
    np.testing.assert_allclose(nt.integrate(ax, ds1.spec, (280, 220)), ref1d["int_metal"], rtol=RTOL)
    np.testing.assert_allclose(nt.mean_ppm(ax, ds1.spec, (280, 220)), ref1d["mean_ppm_metal"], rtol=RTOL)


@pytest.mark.data
def test_integration_with_nmrglue_axis_is_close():
    """Switching to nmrglue's axis changes the Li-metal integral by < 1 %."""
    require(RC1)
    ds = nt.load(data_path(RC1), procno=22)
    ref = golden("rc1_1d_procno22")
    assert ds.integrate((280, 220)) == pytest.approx(float(ref["int_metal"]), rel=1e-2)
    assert ds.mean_ppm((280, 220)) == pytest.approx(float(ref["mean_ppm_metal"]), abs=0.05)


@pytest.mark.data
def test_t1_workflow_matches_legacy():
    """vdlist + integration + mono fit reproduce the legacy T1 analysis."""
    require(f"{PEO}/20")
    ref = golden("peo_t1")
    ds = nt.load(data_path(f"{PEO}/20"))
    tau_ms = nt.read_vdlist(ds.path) * 1000 * nt.vd_factor(ds.path)
    np.testing.assert_allclose(tau_ms, ref["tau"], rtol=1e-12)

    ppm_max = float(ref["ppm_of_max"])
    area = nt.integrate(legacy_ppm_axis(ds.dic), ds.spec, (ppm_max - 10, ppm_max + 10))
    np.testing.assert_allclose(area[: ref["signal"].size], ref["signal"], rtol=RTOL)

    fit = nt.exp_fit(tau_ms, area[: tau_ms.size], "mono", "recovery")
    np.testing.assert_allclose([fit.params[k] for k in ("A", "tau", "y0")],
                               ref["T1mono.params"], rtol=1e-6)
    np.testing.assert_allclose([fit.ci95[k] for k in ("A", "tau", "y0")],
                               ref["T1mono.err"], rtol=1e-4)
    assert fit.r2 == pytest.approx(float(ref["T1mono.R2"]), rel=1e-9)


@pytest.mark.data
def test_trim_empty_rows_matches_legacy():
    """Unacquired rows are dropped as by ``remove_empty=True``."""
    require(IMEC_OCV)
    ref = golden("imec_ocv_pseudo2d")
    ds = nt.load(data_path(IMEC_OCV))
    trimmed = nt.trim_empty_rows(ds)
    np.testing.assert_array_equal(trimmed.spec.shape, ref["trimmed.shape"])
    assert trimmed.ppm_ind.size == trimmed.spec.shape[0]
    assert ds.spec.shape[0] > trimmed.spec.shape[0]


@pytest.mark.data
def test_load_fid_drop_matches_legacy():
    """``group_delay='drop'`` reproduces the legacy FID; time uses 1/SW_h."""
    require(f"{DENDRITES}/71999")
    ref = golden("dendrites_cp_popt")
    fid = nt.load_fid(data_path(f"{DENDRITES}/71999"), group_delay="drop")
    assert_summary_matches(fid.data, ref, "fid999.signal")
    assert fid.time.size == fid.data.shape[-1]
    np.testing.assert_allclose(np.diff(fid.time), 1 / fid.p("SW_h"))


@pytest.mark.data
def test_load_fid_fractional_group_delay():
    """The 500 MHz FID (GRPDLY = 68.014) loads with the digital-filter correction."""
    require(f"{PEO}/27")
    fid = nt.load_fid(data_path(f"{PEO}/27"))
    assert fid.p("GRPDLY") % 1 != 0
    assert fid.data.ndim == 1 and np.all(np.isfinite(fid.data))
    assert fid.time.size == fid.data.size
