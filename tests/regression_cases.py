"""
Regression cases for the legacy ``Datos`` class family.

Each case runs the *current* analysis code (``Datos``, ``DatosProcesados``,
``DatosProcesados2D``, ``DatosProcesadosT1``, ``kpopt_utils`` ...) on a real
dataset, the same way the analysis scripts do, and returns a flat dictionary
of numpy arrays describing the result.

The same functions are used twice:

* ``tests/make_golden.py`` runs them once and stores the output in
  ``tests/golden/<case>.npz`` (the reference, or "golden", results).
* ``tests/test_regression_datos.py`` re-runs them and compares against the
  stored reference, so any refactor that changes a number is caught.

Large arrays are not stored whole. :func:`summarize` keeps their shape,
global sums and an evenly strided sample, which is enough to detect any
real change while keeping the reference files small enough for git.

Adding a case
-------------
Write a function that takes no arguments and returns ``dict[str, array]``,
decorate it with ``@case("dataset/relative/path", ...)`` and run
``python tests/make_golden.py <case_name>``.
"""
import os
from types import SimpleNamespace

import numpy as np

from paths import NMRDATA

#: Registry of all cases: ``name -> (function, list of required paths)``.
CASES = {}

#: Root under which the datasets of the cases live.
DATA_ROOT = NMRDATA


def data_path(relative):
    """
    Build the absolute path of a dataset folder.

    Parameters
    ----------
    relative : str
        Path relative to :data:`DATA_ROOT`, with forward slashes,
        e.g. ``"500/2025-08-05_PEO-PTT-solid-electrolyte_VT/20"``.

    Returns
    -------
    str
        Absolute path, always ending with a separator, because several
        legacy classes still concatenate ``directorio + 'pdata/...'``.
    """
    return os.path.join(DATA_ROOT, *relative.split("/")) + os.sep


def case(*required):
    """
    Register a regression case.

    Parameters
    ----------
    *required : str
        Dataset folders (relative to :data:`DATA_ROOT`) the case needs.
        The test is skipped when any of them is missing, so the suite
        still runs on a machine without the data.

    Returns
    -------
    callable
        Decorator that stores the function in :data:`CASES` under its name.
    """
    def decorator(func):
        CASES[func.__name__] = (func, [data_path(r) for r in required])
        return func
    return decorator


def missing_data(name):
    """
    List the dataset folders of a case that are not present on this machine.

    Parameters
    ----------
    name : str
        Case name (a key of :data:`CASES`).

    Returns
    -------
    list of str
        Missing folders; empty when the case can run.
    """
    return [p for p in CASES[name][1] if not os.path.isdir(p)]


def summarize(out, name, arr, max_points=2048, max_rows=16):
    """
    Add a compact but change-sensitive description of an array to ``out``.

    Stored keys (all prefixed with ``name``):

    ``.shape``
        Shape of the array.
    ``.sum``, ``.abssum``, ``.sqsum``
        Global sums over every element; any change anywhere moves them.
    ``.sample``
        Evenly strided sub-array: at most ``max_rows`` rows (for 2D data)
        and ``max_points`` points along the last axis.

    Complex arrays are stored as two entries, ``name.re`` and ``name.im``.

    Parameters
    ----------
    out : dict
        Dictionary that receives the entries (modified in place).
    name : str
        Prefix for the keys.
    arr : array_like
        Array to describe.
    max_points : int, optional
        Maximum number of samples kept along the last axis.
    max_rows : int, optional
        Maximum number of rows kept for 2D arrays.
    """
    arr = np.asarray(arr)
    if np.iscomplexobj(arr):
        summarize(out, f"{name}.re", arr.real, max_points, max_rows)
        summarize(out, f"{name}.im", arr.imag, max_points, max_rows)
        return
    arr = arr.astype(float)
    out[f"{name}.shape"] = np.array(arr.shape)
    out[f"{name}.sum"] = np.array(arr.sum())
    out[f"{name}.abssum"] = np.array(np.abs(arr).sum())
    out[f"{name}.sqsum"] = np.array((arr ** 2).sum())
    sample = arr
    if sample.ndim >= 1 and sample.shape[-1] > max_points:
        sample = sample[..., ::int(np.ceil(sample.shape[-1] / max_points))]
    if sample.ndim == 2 and sample.shape[0] > max_rows:
        sample = sample[::int(np.ceil(sample.shape[0] / max_rows)), :]
    out[f"{name}.sample"] = sample


def acqus_summary(out, acqus, prefix="acqus"):
    """
    Store the acquisition parameters exposed as attributes by ``Acqus``.

    Parameters
    ----------
    out : dict
        Dictionary that receives the entries (modified in place).
    acqus : Acqus
        Legacy parameter object (``datos.acqus``).
    prefix : str, optional
        Key prefix.
    """
    for attr in ("NS", "RG", "TD", "SW", "SFO1", "BF1", "DW", "D1", "P1", "PL1"):
        if hasattr(acqus, attr):
            out[f"{prefix}.{attr}"] = np.array(getattr(acqus, attr))


# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------
RC1 = "400dnp/2026-09-10_3.2mm_LiMetal_LiF_RC1/265"
PEO = "500/2025-08-05_PEO-PTT-solid-electrolyte_VT"
IMEC_OCV = "500/2025-04-29_IMEC_eLi-sym_LP40-OCV/4"
DENDRITES = "400dnp/2025-06-16_3.2mm_IMECdendrites"


@case(RC1)
def rc1_pseudo2d():
    """
    Pseudo-2D Kpopt series (7Li, DNP) read with ``DatosProcesados2D``.

    Covers the 2D reader, both ppm axes, NS·RG normalisation and
    ``Integrar2D`` on the real, absolute and imaginary parts, before and
    after an in-place ``ppmSelect`` crop (as the Kpopt scripts do).
    """
    from Datos import DatosProcesados2D
    d = DatosProcesados2D(data_path(RC1), p_dir=1)
    out = {}
    acqus_summary(out, d.acqus)
    summarize(out, "ppmAxis", d.espectro.ppmAxis)
    summarize(out, "ppmAxisInd", d.espectro.ppmAxisInd)
    summarize(out, "spec", d.espectro.spec)
    out["int_full"] = d.Integrar2D([300, 10])
    out["int_metal"] = d.Integrar2D([280, 220])
    out["int_metal_abs"] = d.Integrar2D([280, 220], absolute=True)
    out["int_metal_imag"] = d.Integrar2D([280, 220], imaginary=True)
    d.espectro.ppmSelect([300, -50])
    summarize(out, "cropped.ppmAxis", d.espectro.ppmAxis)
    out["cropped.int_full"] = d.Integrar2D([300, 10])
    return out


@case(RC1)
def rc1_1d_procno22():
    """
    1D slice extracted in TopSpin to ``pdata/22``, read with
    ``DatosProcesados``.

    Covers ``p_dir`` selection, the 1D ppm axis, ``Integrar`` on the
    real/absolute/imaginary parts and ``Mean_ppm``.
    """
    from Datos import DatosProcesados
    d = DatosProcesados(data_path(RC1), p_dir=22)
    out = {}
    summarize(out, "ppmAxis", d.espectro.ppmAxis)
    summarize(out, "spec", d.espectro.spec)
    out["int_metal"] = d.Integrar([280, 220])
    out["int_metal_abs"] = d.Integrar([280, 220], absolute=True)
    out["int_metal_imag"] = d.Integrar([280, 220], imaginary=True)
    out["mean_ppm_metal"] = d.Mean_ppm([280, 220])
    out["mean_ppm_full"] = d.Mean_ppm([300, 10])
    return out


@case(RC1)
def rc1_klog_bsms():
    """
    bsms field values read from ``Klog_opt`` by ``kpopt_utils.read_bsms_field``.
    """
    from kpopt_utils import read_bsms_field
    return {"bsms": read_bsms_field(data_path(RC1))}


@case(f"{PEO}/20")
def peo_t1():
    """
    Saturation-recovery T1 (pseudo-2D with ``vdlist``) read with
    ``DatosProcesadosT1``, following ``read_T1data_polymer-electrolyte_KBr.py``.

    Covers the delay list, the vd factor, ``get_T1data`` around the
    maximum of the last spectrum and the mono-exponential ``T1fit``.
    """
    from Datos import DatosProcesadosT1
    d = DatosProcesadosT1(data_path(f"{PEO}/20"))
    out = {}
    acqus_summary(out, d.acqus)
    summarize(out, "ppmAxis", d.espectro.ppmAxis)
    summarize(out, "spec", d.espectro.spec)
    out["factor_vd"] = np.array(d.factor_vd)
    out["tau"] = d.tau
    out["signal_full"] = d.signal
    d.espectro.ppmSelect([100, -100])
    ppm = d.espectro.ppmAxis
    ppm_of_max = ppm[np.argmax(d.espectro.real[-1, :])]
    ppm_range = [ppm_of_max - 10, ppm_of_max + 10]
    out["ppm_of_max"] = np.array(ppm_of_max)
    tau, signal = d.get_T1data(ppm_range)
    out["tau_cut"] = tau
    out["signal"] = signal
    d.T1fit(model="mono")
    out["T1mono.params"] = d.T1params
    out["T1mono.err"] = d.T1params_err
    out["T1mono.R2"] = np.array(d.R2)
    return out


@case(f"{PEO}/27")
def peo_1d_19f():
    """
    Variable-temperature 1D 19F spectrum (``readSpec_PolymerElectrolyte_VT.py``).

    Covers 1D processed data for a nucleus other than 7Li. The raw FID of
    this dataset cannot be read by the legacy code (fractional ``GRPDLY``
    on the 500 MHz console); see ``tests/test_known_bugs.py``.
    """
    from Datos import DatosProcesados
    d = DatosProcesados(data_path(f"{PEO}/27"))
    out = {}
    acqus_summary(out, d.acqus)
    out["gamma"] = np.array(d.gamma)
    summarize(out, "ppmAxis", d.espectro.ppmAxis)
    summarize(out, "spec", d.espectro.spec)
    out["int"] = d.Integrar([-10, -150])
    out["mean_ppm"] = d.Mean_ppm([-10, -150])
    return out


@case(IMEC_OCV)
def imec_ocv_pseudo2d():
    """
    Long in-situ pseudo-2D series (128 MB) as in
    ``read_paeudo2Ddata_integrate_IMEC.py``, plus the ``remove_empty``
    option that trims unacquired rows.
    """
    from Datos import DatosProcesados2D
    path = data_path(IMEC_OCV)
    d = DatosProcesados2D(path)
    out = {}
    acqus_summary(out, d.acqus)
    summarize(out, "ppmAxis", d.espectro.ppmAxis)
    summarize(out, "spec", d.espectro.spec)
    d.espectro.ppmSelect([300, 200])
    out["int"] = d.Integrar2D([300, 200])
    trimmed = DatosProcesados2D(path, remove_empty=True)
    out["trimmed.shape"] = np.array(trimmed.espectro.spec.shape)
    summarize(out, "trimmed.ppmAxisInd", trimmed.espectro.ppmAxisInd)
    return out


@case(f"{DENDRITES}/71", f"{DENDRITES}/71999", f"{DENDRITES}/71998")
def dendrites_cp_popt():
    """
    CP contact-time series stored by popt in ``99n`` sub-experiments, as in
    ``read_paeudo2Ddata_cp-multiple-99n-subfiles.py``.

    Covers ``find_popt_subfolders``, ``extract_popt_parameters``, the
    small 2D reader and the 2D raw-FID (``ser``) reader.
    """
    from Datos import Datos, DatosProcesados2D
    from kpopt_utils import extract_popt_parameters, find_popt_subfolders
    root = data_path(DENDRITES)
    out = {"subfolders": find_popt_subfolders(71, root)}
    for sub in (999, 998):
        out[f"P15_{sub}"] = extract_popt_parameters(sub, "P 15", path=root, expn=71)
        d = DatosProcesados2D(data_path(f"{DENDRITES}/71{sub}"), read_pp=False)
        summarize(out, f"{sub}.ppmAxis", d.espectro.ppmAxis)
        summarize(out, f"{sub}.spec", d.espectro.spec)
        d.espectro.ppmSelect([50, -50])
        out[f"{sub}.int"] = d.Integrar2D([15, -15])
    raw = Datos(data_path(f"{DENDRITES}/71999"), set_fid=True, read_pp=False, ndim=2)
    summarize(out, "fid999.time", raw.fid.timeAxis)
    summarize(out, "fid999.signal", raw.fid.signal)
    return out


@case()
def synthetic_bvalues():
    """
    STE diffusion b-values from ``DatosProcesadosDiff.Crear_bvalue``.

    No diffusion dataset is available on disk, so the method is called on a
    stand-in object with fixed, realistic 7Li parameters (both gradient
    shapes, with and without a calibration factor).
    """
    from Datos import DatosProcesadosDiff
    out = {}
    for shape in ("sin", "rect"):
        for factor_b in (1.0, 1.3):
            fake = SimpleNamespace(
                acqus=SimpleNamespace(gp=80.0),
                acqu2s=SimpleNamespace(TD=16),
                delta=2.0,          # ms
                bigDelta=50.0,      # ms
                gamma=103.962e6,    # 7Li, rad/(s*T)
                factor_b=factor_b,
            )
            DatosProcesadosDiff.Crear_bvalue(fake, shape)
            out[f"{shape}.{factor_b}.bvalue"] = fake.bvalue
            out[f"{shape}.{factor_b}.gpvalue"] = fake.gpvalue
    return out


@case()
def synthetic_gamma():
    """
    Gyromagnetic ratios returned by ``Datos.set_gamma`` for every nucleus
    in its table.
    """
    from Datos import Datos
    nuclei = ["1H", "6Li", "7Li", "13C", "19F", "23Na", "31P", "29Si", "79Br"]
    gammas = [Datos.set_gamma(SimpleNamespace(nucleo=n)) for n in nuclei]
    return {"gamma": np.array(gammas)}
