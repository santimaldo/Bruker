"""
nmrtools: loading and analysing Bruker NMR data on top of nmrglue.

This module is the core of the analysis scripts in this repository. It keeps
the lab's conventions in one place (normalisation by NS·RG, delay lists,
diffusion b-values, gyromagnetic ratios, integration on descending ppm axes)
and otherwise relies on nmrglue's data model:

* :func:`load` returns a :class:`Dataset`: the full nmrglue parameter
  dictionary, the processed spectrum (complex, stored once) and its ppm axes.
* :func:`load_fid` returns a :class:`Fid` with the raw time-domain data.
* Analysis is done with plain functions (:func:`crop`, :func:`integrate`,
  :func:`mean_ppm`, :func:`exp_fit`, ...) that take arrays and return new
  arrays or result objects. Nothing modifies a dataset in place.

Quick start
-----------
>>> import nmrtools as nt
>>> ds = nt.load(r"C:/data/400dnp/2026-09-10_sample/265", procno=1)
>>> ds.p("NS"), ds.p("D1"), ds.pulprog
>>> metal = ds.crop((300, 200))            # new Dataset, ds is untouched
>>> area = metal.integrate((280, 220))     # one value per row for pseudo-2D
>>> fit = nt.exp_fit(nt.read_vdlist(ds.path), area, kind="recovery")
>>> fit.params["tau"], fit.ci95["tau"]

Conventions
-----------
* ppm axes come from nmrglue (``uc_from_udic(...).ppm_scale()``), which
  matches Bruker's definition ``OFFSET - i * SW_p / (SF * SI)``.
* Spectra and FIDs are divided by ``NS * RG`` by default, so experiments
  with different scans or receiver gain can be compared directly.
* Times are in seconds, b-values in s/m^2 and diffusion coefficients in
  m^2/s (SI units), unless a function says otherwise.
* For pseudo-2D data the first axis is the indirect (row) dimension and the
  last axis is the direct (ppm) dimension. Functions that act along ppm
  always use the last axis, so they work the same on 1D and pseudo-2D data.
"""
import os
import re
from dataclasses import dataclass, replace

import nmrglue as ng
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import curve_fit
from scipy.stats import linregress

__all__ = [
    "Dataset", "Fid", "FitResult",
    "load", "load_fid", "trim_empty_rows",
    "crop", "integrate", "mean_ppm",
    "gamma", "GAMMA",
    "read_vdlist", "vd_factor",
    "bvalues", "ste_bvalues",
    "exp_fit", "diffusion_fit",
]

#: Gyromagnetic ratios in rad s^-1 T^-1, keyed by lower-case nucleus label.
#: Source: https://en.wikipedia.org/wiki/Gyromagnetic_ratio
GAMMA = {
    "1h": 267.52218744e6,
    "2h": 41.065e6,
    "6li": 39.994e6,
    "7li": 103.962e6,
    "13c": 67.2828e6,
    "15n": -27.116e6,
    "19f": 251.815e6,
    "23na": 70.761e6,
    "27al": 69.763e6,
    "29si": -53.190e6,
    "31p": 108.291e6,
    "79br": 6.7256e6,
}


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
def _clean_label(value):
    """Strip the ``<...>`` that Bruker puts around string parameters."""
    return str(value).strip().strip("<>").strip()


def gamma(nucleus):
    """
    Gyromagnetic ratio of a nucleus.

    Parameters
    ----------
    nucleus : str
        Nucleus label as Bruker writes it, e.g. ``"7Li"``, ``"<19F>"``,
        ``"1H"``. Case-insensitive.

    Returns
    -------
    float
        Gyromagnetic ratio in rad s^-1 T^-1.

    Raises
    ------
    KeyError
        If the nucleus is not in :data:`GAMMA`. Add it there.
    """
    key = _clean_label(nucleus).lower()
    if key not in GAMMA:
        raise KeyError(f"nucleus '{nucleus}' not in nmrtools.GAMMA; add it there")
    return GAMMA[key]


def _param(dic, name, dim=1):
    """
    Look up a Bruker parameter in an nmrglue dictionary.

    Parameters
    ----------
    dic : dict
        nmrglue dictionary with ``acqus``/``procs`` (and ``acqu2s``/``proc2s``
        for 2D) sub-dictionaries.
    name : str
        Parameter name. Plain names (``"NS"``, ``"SW_h"``, ``"SFO1"``,
        ``"OFFSET"``) are looked up directly. Names ending in a number that
        are not themselves keys (``"D1"``, ``"P1"``, ``"PLW1"``, ``"GPZ6"``,
        ``"CNST2"``) are read as an index into the array parameter, so
        ``"D1"`` is ``dic["acqus"]["D"][1]``.
    dim : int, optional
        1 for the direct dimension (``acqus``, ``procs``), 2 for the indirect
        one (``acqu2s``, ``proc2s``).

    Returns
    -------
    object
        The parameter value (number, string or list).

    Raises
    ------
    KeyError
        If the parameter is not found in the files of that dimension.
    """
    files = {1: ("acqus", "procs"), 2: ("acqu2s", "proc2s")}[dim]
    sources = [dic[f] for f in files if f in dic]
    for src in sources:
        if name in src:
            return src[name]
    match = re.fullmatch(r"([A-Za-z_]+?)(\d+)", name)
    if match:
        base, index = match.group(1), int(match.group(2))
        for src in sources:
            if isinstance(src.get(base), list) and index < len(src[base]):
                return src[base][index]
    raise KeyError(f"parameter '{name}' not found in {files}")


class _ParamMixin:
    """Parameter shortcuts shared by :class:`Dataset` and :class:`Fid`."""

    def p(self, name, dim=1):
        """
        Read a Bruker parameter.

        Parameters
        ----------
        name : str
            Parameter name, e.g. ``"NS"``, ``"RG"``, ``"D1"``, ``"P1"``,
            ``"SW_h"``. See :func:`_param` for how indexed names work.
        dim : int, optional
            1 (default) for the direct dimension, 2 for the indirect one.

        Returns
        -------
        object
            The parameter value.
        """
        return _param(self.dic, name, dim)

    @property
    def nucleus(self):
        """str: Observed nucleus of the direct dimension, e.g. ``"7Li"``."""
        return _clean_label(self.dic["acqus"]["NUC1"])

    @property
    def pulprog(self):
        """str: Name of the pulse program (``PULPROG``)."""
        return _clean_label(self.dic["acqus"]["PULPROG"])

    @property
    def gamma(self):
        """float: Gyromagnetic ratio of :attr:`nucleus`, in rad s^-1 T^-1."""
        return gamma(self.nucleus)


def _norm_factor(dic, normalize):
    """
    Factor that raw intensities are divided by.

    Parameters
    ----------
    dic : dict
        nmrglue dictionary containing ``acqus``.
    normalize : {"ns_rg", "rg", None}
        ``"ns_rg"`` divides by NS·RG, ``"rg"`` only by RG, ``None`` leaves
        the data untouched.

    Returns
    -------
    float
        The divisor.
    """
    if normalize is None:
        return 1.0
    acq = dic["acqus"]
    if normalize == "ns_rg":
        return float(acq["NS"]) * float(acq["RG"])
    if normalize == "rg":
        return float(acq["RG"])
    raise ValueError("normalize must be 'ns_rg', 'rg' or None")


# ---------------------------------------------------------------------------
# Processed data
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Dataset(_ParamMixin):
    """
    A processed Bruker spectrum (1D or pseudo-2D/2D).

    Instances are immutable: methods such as :meth:`crop` return a new
    ``Dataset`` and never modify the original, so the same loaded data can
    be analysed in several ways without reloading.

    Attributes
    ----------
    path : str
        Experiment folder (the one that contains ``acqus``).
    procno : int
        Processing number (``pdata/<procno>``).
    dic : dict
        Full nmrglue parameter dictionary: ``acqus``, ``procs`` and, for 2D
        data, ``acqu2s`` and ``proc2s``. Use :meth:`p` for quick access.
    spec : numpy.ndarray
        Complex spectrum, already divided by :attr:`norm`. Shape ``(n,)``
        for 1D or ``(rows, n)`` for 2D; ppm is always the last axis.
    ppm : numpy.ndarray
        Direct-dimension chemical-shift axis, in ppm (descending).
    ppm_ind : numpy.ndarray or None
        Indirect-dimension axis for 2D data, ``None`` for 1D. For pseudo-2D
        experiments (T1, Kpopt, in-situ series) it carries no physical
        meaning; use the delay list or row index instead.
    title : str
        Contents of ``pdata/<procno>/title`` (empty if there is none).
    norm : float
        Factor the stored intensities were divided by (e.g. NS·RG).
    """

    path: str
    procno: int
    dic: dict
    spec: np.ndarray
    ppm: np.ndarray
    ppm_ind: object
    title: str
    norm: float

    @property
    def real(self):
        """numpy.ndarray: Real part of :attr:`spec`."""
        return self.spec.real

    @property
    def imag(self):
        """numpy.ndarray: Imaginary part of :attr:`spec`."""
        return self.spec.imag

    @property
    def ndim(self):
        """int: Number of dimensions of :attr:`spec` (1 or 2)."""
        return self.spec.ndim

    def crop(self, ppm_range):
        """
        Keep only the points inside a ppm range.

        Parameters
        ----------
        ppm_range : tuple of float
            Two limits in ppm, in any order. Both limits are inclusive.

        Returns
        -------
        Dataset
            New dataset with cropped :attr:`spec` and :attr:`ppm`.
        """
        ppm, spec = crop(self.ppm, self.spec, ppm_range)
        return replace(self, ppm=ppm, spec=spec)

    def rows(self, selection):
        """
        Select rows of a 2D dataset (the indirect dimension).

        Parameters
        ----------
        selection : slice, int array or bool array
            Any numpy index for the first axis, e.g. ``slice(1, None)`` to
            drop the first row.

        Returns
        -------
        Dataset
            New dataset with the selected rows (and matching ``ppm_ind``).
        """
        if self.ndim != 2:
            raise ValueError("rows() needs a 2D dataset")
        spec = self.spec[selection, :]
        ppm_ind = None if self.ppm_ind is None else self.ppm_ind[selection]
        return replace(self, spec=spec, ppm_ind=ppm_ind)

    def integrate(self, ppm_range=None, part="real"):
        """
        Integrate the spectrum over a ppm range. See :func:`integrate`.

        Parameters
        ----------
        ppm_range : tuple of float, optional
            Integration limits; whole axis if omitted.
        part : {"real", "imag", "abs", "complex"}, optional
            Which part of the spectrum to integrate.

        Returns
        -------
        float or numpy.ndarray
            One value for 1D data, one value per row for 2D data.
        """
        return integrate(self.ppm, self.spec, ppm_range, part)

    def mean_ppm(self, ppm_range=None, part="real"):
        """
        Intensity-weighted mean chemical shift. See :func:`mean_ppm`.

        Parameters
        ----------
        ppm_range : tuple of float, optional
            Limits of the region; whole axis if omitted.
        part : {"real", "imag", "abs"}, optional
            Which part of the spectrum is used as weight.

        Returns
        -------
        float or numpy.ndarray
            One value for 1D data, one value per row for 2D data.
        """
        return mean_ppm(self.ppm, self.spec, ppm_range, part)


def _pdata_files(pdir):
    """
    Names of the real and imaginary binary files present in a pdata folder.

    Parameters
    ----------
    pdir : str
        ``<expno>/pdata/<procno>`` folder.

    Returns
    -------
    tuple of (str, str or None)
        Real file (``1r``/``2rr``) and imaginary file (``1i``/``2ii``, or
        ``None`` when it was not stored).
    """
    for real, imag in (("1r", "1i"), ("2rr", "2ii")):
        if os.path.isfile(os.path.join(pdir, real)):
            return real, (imag if os.path.isfile(os.path.join(pdir, imag)) else None)
    raise FileNotFoundError(f"no 1r or 2rr file in {pdir}")


def _ppm_axis(dic, data, axis):
    """
    Chemical-shift axis of one dimension, computed by nmrglue.

    Parameters
    ----------
    dic : dict
        nmrglue dictionary with acqus/procs (and acqu2s/proc2s).
    data : numpy.ndarray
        Real data array, used for the dimension sizes.
    axis : int
        numpy axis of ``data`` (``-1`` direct, ``0`` indirect for 2D).

    Returns
    -------
    numpy.ndarray
        ppm values, one per point along ``axis``.
    """
    procs = dic["procs"]
    stsi = procs.get("STSI", 0)
    strip_fake = bool(stsi) and stsi != procs.get("SI")
    udic = ng.bruker.guess_udic(dic, data, strip_fake=strip_fake)
    return ng.fileiobase.uc_from_udic(udic, axis % data.ndim).ppm_scale()


def load(path, procno=1, normalize="ns_rg"):
    """
    Load a processed Bruker spectrum (``1r/1i`` or ``2rr/2ii``).

    Parameters
    ----------
    path : str
        Experiment folder, i.e. the folder that contains ``acqus`` and
        ``pdata``. A trailing separator is optional.
    procno : int, optional
        Processing number to read from ``pdata/<procno>``. Default 1.
    normalize : {"ns_rg", "rg", None}, optional
        Divide intensities by NS·RG (default), by RG only, or not at all.

    Returns
    -------
    Dataset
        The spectrum with its parameters and ppm axes.

    Notes
    -----
    The real and imaginary files are requested by name, so a ``2ri`` or
    ``2ir`` file in the folder can never be mistaken for the imaginary part.
    If no imaginary file was stored, the imaginary part is zero.
    """
    path = os.path.normpath(path)
    pdir = os.path.join(path, "pdata", str(procno))
    real_file, imag_file = _pdata_files(pdir)
    files = [real_file] + ([imag_file] if imag_file else [])
    dic, parts = ng.bruker.read_pdata(pdir, bin_files=files, all_components=True)
    real = np.asarray(parts[0], dtype=float)
    imag = np.asarray(parts[1], dtype=float) if imag_file else np.zeros_like(real)

    norm = _norm_factor(dic, normalize)
    spec = (real + 1j * imag) / norm
    ppm = _ppm_axis(dic, real, -1)
    ppm_ind = _ppm_axis(dic, real, 0) if real.ndim == 2 else None

    title_file = os.path.join(pdir, "title")
    title = ""
    if os.path.isfile(title_file):
        with open(title_file, "r", errors="replace") as f:
            title = f.read().strip()
    return Dataset(path=path, procno=procno, dic=dic, spec=spec, ppm=ppm,
                   ppm_ind=ppm_ind, title=title, norm=norm)


def trim_empty_rows(ds, threshold=0.0):
    """
    Drop trailing rows of a 2D dataset that were never acquired.

    TopSpin pads a pseudo-2D series to the processed size (SI of F1) even if
    the experiment was stopped early; those rows are all zero.

    Parameters
    ----------
    ds : Dataset
        2D dataset.
    threshold : float, optional
        A row counts as empty when every ``|spec|`` value is ``<= threshold``.

    Returns
    -------
    Dataset
        New dataset ending at the last non-empty row. Unchanged if all rows
        contain data (or if all are empty).
    """
    if ds.ndim != 2:
        raise ValueError("trim_empty_rows() needs a 2D dataset")
    filled = np.any(np.abs(ds.spec) > threshold, axis=1)
    if not filled.any():
        return ds
    last = int(np.nonzero(filled)[0][-1]) + 1
    return ds.rows(slice(0, last))


# ---------------------------------------------------------------------------
# Raw data
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Fid(_ParamMixin):
    """
    Raw Bruker time-domain data (``fid`` or ``ser``).

    Attributes
    ----------
    path : str
        Experiment folder.
    dic : dict
        Full nmrglue parameter dictionary (``acqus``, and ``acqu2s`` for 2D).
    data : numpy.ndarray
        Complex FID(s), divided by :attr:`norm`. Shape ``(n,)`` or
        ``(rows, n)``; time is always the last axis.
    time : numpy.ndarray
        Time of each complex point, in seconds, starting at 0. The spacing
        is ``1 / SW_h`` (one complex point every two Bruker ``DW``).
    norm : float
        Factor the stored intensities were divided by (e.g. NS·RG).
    group_delay : str or None
        How the digital-filter group delay was handled (see :func:`load_fid`).
    """

    path: str
    dic: dict
    data: np.ndarray
    time: np.ndarray
    norm: float
    group_delay: object


def load_fid(path, normalize="ns_rg", group_delay="filter"):
    """
    Load raw Bruker time-domain data.

    Parameters
    ----------
    path : str
        Experiment folder containing ``fid`` or ``ser``.
    normalize : {"ns_rg", "rg", None}, optional
        Divide intensities by NS·RG (default), by RG only, or not at all.
    group_delay : {"filter", "drop", None}, optional
        How to deal with the points delayed by Bruker's digital filter:

        ``"filter"`` (default)
            nmrglue's :func:`nmrglue.bruker.remove_digital_filter`: a
            frequency-domain first-order phase shift by the exact (possibly
            fractional) ``GRPDLY``. This is what TopSpin does.
        ``"drop"``
            Discard the first ``int(GRPDLY) + 1`` points, as the legacy
            ``Datos.set_fid`` did. Kept for reproducing old results; it also
            works on consoles with a fractional ``GRPDLY``.
        ``None``
            Leave the data untouched.

    Returns
    -------
    Fid
        The raw data with its time axis and parameters.

    Notes
    -----
    nmrglue cannot always infer the row layout of a pseudo-2D ``ser`` file
    and then returns 1D data. In that case the data is split into rows of
    ``TD/2`` complex points, discarding an incomplete last row.
    """
    path = os.path.normpath(path)
    dic, data = ng.bruker.read(path)
    data = np.asarray(data)
    npts = int(dic["acqus"]["TD"]) // 2
    if "acqu2s" in dic and data.ndim == 1 and data.size > npts:
        rows = data.size // npts
        data = data[: rows * npts].reshape(rows, npts)

    if group_delay == "filter":
        data = ng.bruker.remove_digital_filter(dic, data, truncate=False)
    elif group_delay == "drop":
        data = data[..., int(dic["acqus"]["GRPDLY"]) + 1:]
    elif group_delay is not None:
        raise ValueError("group_delay must be 'filter', 'drop' or None")

    norm = _norm_factor(dic, normalize)
    data = data / norm
    time = np.arange(data.shape[-1]) / float(dic["acqus"]["SW_h"])
    return Fid(path=path, dic=dic, data=data, time=time, norm=norm,
               group_delay=group_delay)


# ---------------------------------------------------------------------------
# Spectral analysis (pure functions on arrays)
# ---------------------------------------------------------------------------
def _part(spec, part):
    """
    Select the real, imaginary, absolute or complex part of a spectrum.

    Parameters
    ----------
    spec : numpy.ndarray
        Complex (or real) spectrum.
    part : {"real", "imag", "abs", "complex"}
        Part to return.

    Returns
    -------
    numpy.ndarray
    """
    if part == "real":
        return np.real(spec)
    if part == "imag":
        return np.imag(spec)
    if part == "abs":
        return np.abs(spec)
    if part == "complex":
        return np.asarray(spec)
    raise ValueError("part must be 'real', 'imag', 'abs' or 'complex'")


def _inside(ppm, ppm_range):
    """
    Boolean mask of the points strictly inside a ppm range.

    Parameters
    ----------
    ppm : numpy.ndarray
        ppm axis.
    ppm_range : tuple of float or None
        Two limits in any order; ``None`` selects every point.

    Returns
    -------
    numpy.ndarray of bool
    """
    if ppm_range is None:
        return np.ones(ppm.shape, dtype=bool)
    lo, hi = min(ppm_range), max(ppm_range)
    return np.abs(ppm - (lo + hi) / 2) < (hi - lo) / 2


def crop(ppm, spec, ppm_range):
    """
    Cut a spectrum (1D or 2D) to a ppm range.

    Parameters
    ----------
    ppm : numpy.ndarray
        ppm axis, same length as the last axis of ``spec``.
    spec : numpy.ndarray
        Spectrum, shape ``(n,)`` or ``(rows, n)``.
    ppm_range : tuple of float
        Two limits in any order. Both limits are inclusive.

    Returns
    -------
    ppm : numpy.ndarray
        Cropped axis.
    spec : numpy.ndarray
        Cropped spectrum (a new array; the input is not modified).

    Raises
    ------
    ValueError
        If no point falls inside the range.
    """
    lo, hi = min(ppm_range), max(ppm_range)
    keep = (ppm >= lo) & (ppm <= hi)
    if not keep.any():
        raise ValueError(f"no points between {lo} and {hi} ppm")
    return ppm[keep].copy(), np.asarray(spec)[..., keep].copy()


def integrate(ppm, spec, ppm_range=None, part="real"):
    """
    Integrate a spectrum over a ppm range with Simpson's rule.

    The sign is independent of the axis direction: a positive peak gives a
    positive integral on Bruker's descending ppm axis.

    Parameters
    ----------
    ppm : numpy.ndarray
        ppm axis, same length as the last axis of ``spec``.
    spec : numpy.ndarray
        Spectrum, shape ``(n,)`` or ``(rows, n)``.
    ppm_range : tuple of float, optional
        Integration limits in any order. Only points strictly inside the
        limits are used (as in the legacy ``Datos.Integrar``). Whole axis if
        omitted.
    part : {"real", "imag", "abs", "complex"}, optional
        Part of the spectrum to integrate. Default ``"real"``.

    Returns
    -------
    float or numpy.ndarray
        Integral in (intensity · ppm). One value for 1D input, one per row
        for 2D input.
    """
    mask = _inside(ppm, ppm_range)
    x = ppm[mask]
    y = _part(spec, part)[..., mask]
    area = simpson(y, x=x, axis=-1)
    return -area if x[0] > x[-1] else area


def mean_ppm(ppm, spec, ppm_range=None, part="real"):
    """
    Intensity-weighted mean chemical shift (centre of mass) of a region.

    Computes ``∫ ppm·S(ppm) dppm / ∫ S(ppm) dppm`` over the region.

    Parameters
    ----------
    ppm : numpy.ndarray
        ppm axis, same length as the last axis of ``spec``.
    spec : numpy.ndarray
        Spectrum, shape ``(n,)`` or ``(rows, n)``.
    ppm_range : tuple of float, optional
        Limits of the region (points strictly inside are used). Whole axis
        if omitted.
    part : {"real", "imag", "abs"}, optional
        Part of the spectrum used as weight. Default ``"real"``.

    Returns
    -------
    float or numpy.ndarray
        Mean position in ppm; one value per row for 2D input.
    """
    mask = _inside(ppm, ppm_range)
    x = ppm[mask]
    y = _part(spec, part)[..., mask]
    return simpson(x * y, x=x, axis=-1) / simpson(y, x=x, axis=-1)


# ---------------------------------------------------------------------------
# Delay lists and diffusion
# ---------------------------------------------------------------------------
_TIME_UNITS = {"u": 1e-6, "m": 1e-3, "s": 1.0}


def read_vdlist(path, filename="vdlist"):
    """
    Read a TopSpin variable-delay list.

    Parameters
    ----------
    path : str
        Experiment folder (or any folder) containing the list.
    filename : str, optional
        Name of the list file, relative to ``path``. Default ``"vdlist"``.

    Returns
    -------
    numpy.ndarray
        Delays in seconds. Values may carry TopSpin's unit suffixes
        (``u``, ``m``, ``s``); plain numbers are taken as seconds.
    """
    values = []
    with open(os.path.join(path, filename), "r") as f:
        for line in f:
            token = line.split("#")[0].strip()
            if not token:
                continue
            unit = token[-1].lower()
            if unit in _TIME_UNITS:
                values.append(float(token[:-1]) * _TIME_UNITS[unit])
            else:
                values.append(float(token))
    return np.array(values)


def vd_factor(path, default=1.0):
    """
    Multiplier applied to ``vd`` inside the pulse program (``vd*N``).

    Some sequences use the listed delay several times (e.g. ``vd*2``), so the
    real delay is the vdlist value times that factor.

    Parameters
    ----------
    path : str
        Experiment folder containing ``pulseprogram``.
    default : float, optional
        Value returned when no ``vd*N`` line is found. Default 1.

    Returns
    -------
    float
        The factor ``N``.
    """
    with open(os.path.join(path, "pulseprogram"), "r", errors="replace") as f:
        for line in f:
            match = re.match(r"\s*vd\s*\*\s*([0-9.eE+-]+)", line)
            if match:
                return float(match.group(1))
    return default


def bvalues(gp_percent, delta, big_delta, gamma, g_max=12.0, shape="sin"):
    """
    Diffusion weighting b for a stimulated-echo (STE) experiment.

    Parameters
    ----------
    gp_percent : array_like
        Gradient amplitudes in % of ``g_max``.
    delta : float
        Gradient pulse length δ in s (total, e.g. ``2*P30`` for bipolar).
    big_delta : float
        Diffusion time Δ in s.
    gamma : float
        Gyromagnetic ratio in rad s^-1 T^-1 (see :func:`gamma`).
    g_max : float, optional
        Gradient strength at 100 %, in T/m. Default 12.
    shape : {"sin", "rect"}, optional
        Gradient pulse shape:

        * ``"sin"``: ``b = (γ g δ / π)^2 (4Δ - δ)``
        * ``"rect"``: ``b = (γ g δ)^2 (Δ - δ/3)``

    Returns
    -------
    numpy.ndarray
        b-values in s/m^2.
    """
    g = np.asarray(gp_percent, dtype=float) / 100.0 * g_max
    if shape == "sin":
        return (gamma * g * delta / np.pi) ** 2 * (4 * big_delta - delta)
    if shape == "rect":
        return (gamma * g * delta) ** 2 * (big_delta - delta / 3)
    raise ValueError("shape must be 'sin' or 'rect'")


def ste_bvalues(data, g_max=12.0, shape="sin"):
    """
    Gradient list and b-values of a bipolar STE diffusion pseudo-2D.

    Reads the sequence parameters the way the legacy ``DatosProcesadosDiff``
    did: δ = 2·P30, Δ = D20, maximum gradient GPZ6, and a linear ramp from
    0 to GPZ6 over the TD points of the indirect dimension.

    Parameters
    ----------
    data : Dataset or Fid
        Any object with a ``p()`` parameter method and a 2D ``dic``.
    g_max : float, optional
        Gradient strength at 100 %, in T/m. Default 12.
    shape : {"sin", "rect"}, optional
        Gradient pulse shape, see :func:`bvalues`.

    Returns
    -------
    gp : numpy.ndarray
        Gradient amplitudes in %.
    b : numpy.ndarray
        b-values in s/m^2.
    """
    delta = 2 * data.p("P30") * 1e-6
    big_delta = data.p("D20")
    gp = np.linspace(0, data.p("GPZ6"), int(data.p("TD", dim=2)))
    return gp, bvalues(gp, delta, big_delta, data.gamma, g_max, shape)


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FitResult:
    """
    Result of a least-squares fit.

    Attributes
    ----------
    model : str
        Name of the model that was fitted.
    params : dict of str to float
        Best-fit parameters.
    stderr : dict of str to float
        One-standard-deviation uncertainties.
    r2 : float
        Coefficient of determination of the fit, on the fitted data.
    x : numpy.ndarray
        x data that was fitted (sorted ascending).
    y : numpy.ndarray
        y data that was fitted (in the same order as :attr:`x`).
    func : callable
        Model function ``func(x, *values)``, in the order of :attr:`params`.
    """

    model: str
    params: dict
    stderr: dict
    r2: float
    x: np.ndarray
    y: np.ndarray
    func: object

    @property
    def ci95(self):
        """dict of str to float: 95 % confidence half-widths (1.96·stderr)."""
        return {k: 1.96 * v for k, v in self.stderr.items()}

    @property
    def residuals(self):
        """numpy.ndarray: ``y - predict(x)`` on the fitted data."""
        return self.y - self.predict(self.x)

    def predict(self, x):
        """
        Evaluate the fitted model.

        Parameters
        ----------
        x : array_like
            Points at which to evaluate.

        Returns
        -------
        numpy.ndarray
        """
        return self.func(np.asarray(x, dtype=float), *self.params.values())

    def curve(self, n=500, log=None):
        """
        Smooth curve of the fitted model over the fitted x range, for plots.

        Parameters
        ----------
        n : int, optional
            Number of points. Default 500.
        log : bool, optional
            Logarithmic spacing. By default it is used when all x > 0 and
            they span more than two decades.

        Returns
        -------
        x : numpy.ndarray
        y : numpy.ndarray
        """
        lo, hi = self.x.min(), self.x.max()
        if log is None:
            log = lo > 0 and hi / lo > 100
        x = np.geomspace(lo, hi, n) if log else np.linspace(lo, hi, n)
        return x, self.predict(x)


def _r2(y, y_fit):
    """Coefficient of determination."""
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1 - ss_res / ss_tot)


def _mono(x, A, tau, y0):
    """Mono-exponential model ``y0 + A exp(-x/tau)``."""
    return y0 + A * np.exp(-x / tau)


def _bi(x, A1, tau1, A2, tau2, y0):
    """Bi-exponential model ``y0 + A1 exp(-x/tau1) + A2 exp(-x/tau2)``."""
    return y0 + A1 * np.exp(-x / tau1) + A2 * np.exp(-x / tau2)


def exp_fit(x, y, model="mono", kind="recovery", p0=None, bounds=None):
    """
    Fit a mono- or bi-exponential, e.g. T1 recovery or T2 decay.

    Models:

    * ``"mono"``: ``y0 + A exp(-x/tau)``, parameters ``A, tau, y0``.
    * ``"bi"``: ``y0 + A1 exp(-x/tau1) + A2 exp(-x/tau2)``, parameters
      ``A1, tau1, A2, tau2, y0``.

    Parameters
    ----------
    x : array_like
        Delays (e.g. from :func:`read_vdlist`), any units; ``tau`` comes out
        in the same units.
    y : array_like
        Signal, e.g. integrals from :func:`integrate`.
    model : {"mono", "bi"}, optional
        Number of exponential components. Default ``"mono"``.
    kind : {"recovery", "decay"}, optional
        Constrains the sign of the amplitudes: ``"recovery"`` (saturation
        or inversion recovery, signal grows with x) forces ``A <= 0``;
        ``"decay"`` (T2, CPMG) forces ``A >= 0``. Default ``"recovery"``.
    p0 : sequence of float, optional
        Initial guess in parameter order. Estimated from the data if omitted.
    bounds : tuple of sequences, optional
        ``(lower, upper)`` bounds in parameter order, overriding ``kind``.

    Returns
    -------
    FitResult
        Parameters, standard errors, R², and the model for plotting.

    Raises
    ------
    RuntimeError
        If the least-squares fit does not converge.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    if kind not in ("recovery", "decay"):
        raise ValueError("kind must be 'recovery' or 'decay'")
    sign = -1.0 if kind == "recovery" else 1.0

    # Rough guess: plateau = last point, amplitude = first - last, and the
    # time constant where the signal has covered 1 - 1/e of the change.
    y0_guess = y[-1]
    a_guess = y[0] - y[-1]
    target = y0_guess + a_guess / np.e
    tau_guess = float(x[np.argmin(np.abs(y - target))]) or float(np.median(x[x > 0]))

    if model == "mono":
        func, names = _mono, ["A", "tau", "y0"]
        guess = [a_guess, tau_guess, y0_guess]
        amp = [0]
    elif model == "bi":
        func, names = _bi, ["A1", "tau1", "A2", "tau2", "y0"]
        guess = [0.5 * a_guess, tau_guess / 3, 0.5 * a_guess, tau_guess * 3, y0_guess]
        amp = [0, 2]
    else:
        raise ValueError("model must be 'mono' or 'bi'")

    if bounds is None:
        lower = [-np.inf] * len(names)
        upper = [np.inf] * len(names)
        for i in amp:
            if sign < 0:
                upper[i] = 0.0
            else:
                lower[i] = 0.0
        for i, name in enumerate(names):
            if name.startswith("tau"):
                lower[i] = 0.0
        bounds = (lower, upper)
    if p0 is None:
        # A guess on a bound is fine: scipy moves it strictly inside.
        p0 = np.clip(guess, bounds[0], bounds[1])

    popt, pcov = curve_fit(func, x, y, p0=p0, bounds=bounds)
    err = np.sqrt(np.diag(pcov))
    return FitResult(model=f"exp-{model}-{kind}",
                     params={n: float(v) for n, v in zip(names, popt)},
                     stderr={n: float(e) for n, e in zip(names, err)},
                     r2=_r2(y, func(x, *popt)),
                     x=x, y=y, func=func)


def _diffusion_model(b, D, S0):
    """Stejskal–Tanner attenuation ``S0 exp(-b D)``."""
    return S0 * np.exp(-b * D)


def diffusion_fit(b, signal):
    """
    Diffusion coefficient from a linear fit of ``ln(S)`` against b.

    Parameters
    ----------
    b : array_like
        b-values in s/m^2 (see :func:`bvalues`).
    signal : array_like
        Positive signal intensities at each b.

    Returns
    -------
    FitResult
        ``params["D"]`` in m^2/s and ``params["S0"]``, with standard errors.
        ``r2`` is that of the linear fit of ``ln(S)``.
    """
    b = np.asarray(b, dtype=float)
    s = np.asarray(signal, dtype=float)
    order = np.argsort(b)
    b, s = b[order], s[order]
    reg = linregress(b, np.log(s))
    s0 = float(np.exp(reg.intercept))
    return FitResult(model="diffusion-linear",
                     params={"D": float(-reg.slope), "S0": s0},
                     stderr={"D": float(reg.stderr), "S0": float(s0 * reg.intercept_stderr)},
                     r2=float(reg.rvalue ** 2), x=b, y=s, func=_diffusion_model)
