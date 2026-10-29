"""
Tests for ``relaxation_fit.py`` (the reusable T1/T2 script), for
:mod:`provenance` and for ``nmrtools`` start times.

The script's analysis functions are imported without running it (its run
block is under ``if __name__ == "__main__"``). The data test reproduces the
legacy T1 analysis of ``read_T1data_polymer-electrolyte_KBr.py`` stored in
``tests/golden/peo_t1.npz``.
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import nmrtools as nt  # noqa: E402
import relaxation_fit as rf  # noqa: E402
from provenance import RECORD_NAME, write_run_record  # noqa: E402
from regression_cases import PEO, data_path  # noqa: E402

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------
def test_run_record_contents(tmp_path):
    """The record stores parameters (JSON-converted), script name and git info."""
    params = {"range": (300, 200), "arr": np.arange(3), "x": np.float64(1.5)}
    record = write_run_record(str(tmp_path / "out"), params, rf.__file__)
    with open(tmp_path / "out" / RECORD_NAME, encoding="utf-8") as f:
        saved = json.load(f)
    assert saved == record
    assert saved["script"] == "relaxation_fit.py"
    assert saved["parameters"] == {"range": [300, 200], "arr": [0, 1, 2], "x": 1.5}
    assert set(saved["git"]) == {"commit", "dirty"}


def test_run_record_snapshot_when_dirty(tmp_path):
    """A script outside any git repository is copied next to the record."""
    script = tmp_path / "analysis.py"
    script.write_text("print('hi')\n")
    record = write_run_record(str(tmp_path / "out"), {}, str(script))
    assert record["git"]["dirty"] is True
    assert (tmp_path / "out" / "analysis.py.snapshot").read_text() == "print('hi')\n"


# ---------------------------------------------------------------------------
# helpers of the script
# ---------------------------------------------------------------------------
PPM = np.linspace(100, -100, 2001)


def test_signal_windows_centred_and_fixed():
    """Centred windows follow the maximum of the reference row; fixed ones pass through."""
    spec = np.zeros((3, PPM.size))
    spec[0, np.argmin(np.abs(PPM - 40))] = 1.0      # first row peak at 40 ppm
    spec[-1, np.argmin(np.abs(PPM + 20))] = 1.0     # last row peak at -20 ppm
    win, pmax, _ = rf.signal_windows(PPM, spec, ("centred", 10), kind="recovery")
    assert pmax == pytest.approx(-20) and win == [(pytest.approx(-15), pytest.approx(-25))]
    win, pmax, _ = rf.signal_windows(PPM, spec, ("centred", 10), kind="decay")
    assert pmax == pytest.approx(40)
    win, _, _ = rf.signal_windows(PPM, spec, ("fixed", [(30, 10), (0, -5)]))
    assert win == [(30, 10), (0, -5)]
    with pytest.raises(ValueError):
        rf.signal_windows(PPM, spec, ("around", 5))


def test_fit_row_bi_sorts_components():
    """Bi-exponential rows list the short component first, with fractions."""
    tau = np.geomspace(1e-4, 10, 30)
    y = 100 - 30 * np.exp(-tau / 2.0) - 70 * np.exp(-tau / 0.01)
    row = rf.fit_row(nt.exp_fit(tau, y, "bi", "recovery"), "bi")
    assert row["T1_short_s"] == pytest.approx(0.01, rel=1e-4)
    assert row["T1_long_s"] == pytest.approx(2.0, rel=1e-4)
    assert row["frac_short"] == pytest.approx(0.7, rel=1e-4)
    assert np.isnan(rf.fit_row(None, "bi")["T1_short_s"])


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
@pytest.mark.data
def test_start_time_from_acqus():
    """``start_time`` comes from DATE_START (here 2025-08-05)."""
    if not os.path.isdir(data_path(f"{PEO}/20")):
        pytest.skip("dataset not on this machine")
    ds = nt.load(data_path(f"{PEO}/20"))
    assert ds.start_time.timestamp() == pytest.approx(ds.p("DATE_START"))
    assert ds.start_time.date().isoformat() == "2025-08-05"


@pytest.mark.data
def test_t1_matches_legacy_analysis(monkeypatch):
    """The centred-window mono fit reproduces the legacy KBr-script T1."""
    if not os.path.isdir(data_path(f"{PEO}/20")):
        pytest.skip("dataset not on this machine")
    monkeypatch.setattr(rf, "PATH", data_path(PEO))
    result = rf.analyse_experiment(dict(label="T-20C", expno=20, signal=("centred", 20)))
    with np.load(os.path.join(GOLDEN, "peo_t1.npz")) as ref:
        assert result["ppm_of_max"] == pytest.approx(float(ref["ppm_of_max"]), abs=1e-9)
        np.testing.assert_allclose(result["tau"] * 1000, ref["tau_cut"], rtol=1e-12)
        np.testing.assert_allclose(result["signals"][0], ref["signal"], rtol=1e-9)
        row = rf.fit_row(result["fits"][0], "mono")
        assert row["T1_s"] * 1000 == pytest.approx(float(ref["T1mono.params"][1]), rel=1e-6)
        assert row["T1_ci95_s"] * 1000 == pytest.approx(float(ref["T1mono.err"][1]), rel=1e-4)
        assert row["r2"] == pytest.approx(float(ref["T1mono.R2"]), rel=1e-9)
