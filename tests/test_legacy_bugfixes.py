"""
Regression tests for bugs of the legacy ``Datos`` classes fixed in Phase 2.

Until Phase 2 these tests were marked ``xfail(strict=True)`` and documented
the bugs. Phase 2 made the legacy classes delegate to :mod:`nmrtools`, which
fixed them; the tests now guard against the bugs coming back.
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
def test_integrate_absolute_after_ppmselect_1d():
    """
    ``Integrar(absolute=True)`` works after an in-place ``ppmSelect`` crop.

    Was: ``ppmSelect`` cropped real/imag/ppmAxis but not ``spec`` (1D), so
    the absolute integral raised ``IndexError``. ``Espectro.spec`` is now
    computed from real and imag.
    """
    require(RC1)
    from Datos import DatosProcesados
    d = DatosProcesados(data_path(RC1), p_dir=22)
    expected = d.Integrar([280, 220], absolute=True)
    d.espectro.ppmSelect([300, 200])
    assert d.espectro.spec.shape == d.espectro.ppmAxis.shape
    np.testing.assert_allclose(d.Integrar([280, 220], absolute=True), expected, rtol=1e-6)


@pytest.mark.data
def test_fid_time_axis_spacing():
    """
    ``Datos.set_fid`` spaces complex points by ``1/SW_h``.

    Was: the spacing was Bruker's DW (dwell per *real* point), so the time
    axis was half the real acquisition time.
    """
    require(f"{DENDRITES}/71999")
    from Datos import Datos
    raw = Datos(data_path(f"{DENDRITES}/71999"), set_fid=True, read_pp=False, ndim=2)
    step = raw.fid.timeAxis[1] - raw.fid.timeAxis[0]
    np.testing.assert_allclose(step, 1 / raw.acqus.dic["SW_h"], rtol=1e-9)


@pytest.mark.data
def test_read_fid_with_fractional_group_delay():
    """
    ``Datos(set_fid=True)`` reads a 1D FID from the 500 MHz spectrometer.

    Was: ``GRPDLY`` (68.014 on that console) was used as a slice index and
    raised ``TypeError``. The digital filter is now removed by nmrglue.
    """
    require(f"{PEO}/27")
    from Datos import Datos
    raw = Datos(data_path(f"{PEO}/27"), set_fid=True)
    assert raw.fid.signal.ndim == 1 and raw.fid.signal.size > 0


@pytest.mark.data
def test_pseudo2d_indirect_axis_not_empty():
    """
    ``DatosProcesados2D`` gives one indirect-axis value per row.

    Was: ``proc2s`` of a pseudo-2D has ``FTSIZE = 0``, so ``ppmAxisInd`` and
    the ``ppmGrid*`` meshgrids were empty arrays.
    """
    require(RC1)
    from Datos import DatosProcesados2D
    d = DatosProcesados2D(data_path(RC1), p_dir=1)
    rows, cols = d.espectro.spec.shape
    assert d.espectro.ppmAxisInd.size == rows
    assert d.espectro.ppmGridDir.shape == (rows, cols)
