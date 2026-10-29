"""
Known bugs in the legacy ``Datos`` classes, written as tests of the
*correct* behaviour.

Each test is marked ``xfail(strict=True)``: it is expected to fail today.
When a later refactoring phase fixes the bug, the test starts passing,
pytest reports it as ``XPASS(strict)`` (a failure), and the marker must be
removed. That way a fix can never go unnoticed, and the test then guards
against the bug coming back.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from regression_cases import DENDRITES, PEO, RC1, data_path  # noqa: E402


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


@pytest.mark.data
@pytest.mark.xfail(strict=True, raises=IndexError,
                   reason="ppmSelect crops real/imag/ppmAxis but not spec (1D); "
                          "Integrar(absolute=True) then indexes the uncropped spec")
def test_integrate_absolute_after_ppmselect_1d():
    """``Integrar(absolute=True)`` works after an in-place ``ppmSelect`` crop."""
    require(RC1)
    from Datos import DatosProcesados
    d = DatosProcesados(data_path(RC1), p_dir=22)
    expected = d.Integrar([280, 220], absolute=True)
    d.espectro.ppmSelect([300, 200])
    np.testing.assert_allclose(d.Integrar([280, 220], absolute=True), expected, rtol=1e-6)


@pytest.mark.data
@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="Acqus.DW = 1/(2*SW_h) is the dwell per real point, but set_fid "
                          "uses it as the spacing of complex points: time axis is halved")
def test_fid_time_axis_spacing():
    """``Datos.set_fid`` spaces complex points by ``1/SW_h`` (true acquisition time)."""
    require(f"{DENDRITES}/71999")
    from Datos import Datos
    raw = Datos(data_path(f"{DENDRITES}/71999"), set_fid=True, read_pp=False, ndim=2)
    step = raw.fid.timeAxis[1] - raw.fid.timeAxis[0]
    np.testing.assert_allclose(step, 1 / raw.acqus.dic["SW_h"], rtol=1e-9)


@pytest.mark.data
@pytest.mark.xfail(strict=True, raises=TypeError,
                   reason="set_fid uses GRPDLY as a slice index; the 500 MHz console "
                          "stores a fractional group delay (68.014)")
def test_read_fid_with_fractional_group_delay():
    """``Datos(set_fid=True)`` reads a 1D FID from the 500 MHz spectrometer."""
    require(f"{PEO}/27")
    from Datos import Datos
    raw = Datos(data_path(f"{PEO}/27"), set_fid=True)
    assert raw.fid.signal.ndim == 1 and raw.fid.signal.size > 0
