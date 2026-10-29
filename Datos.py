# coding=UTF-8
"""
Legacy interface to Bruker data: ``Datos`` and its subclasses.

These classes are kept so that the existing analysis scripts keep working
unchanged. Since the nmrtools refactor (Phase 2) they no longer read or
process data themselves: loading, ppm axes, integration, delay lists,
b-values and fits are all delegated to :mod:`nmrtools`. New code should use
:mod:`nmrtools` directly.

Behaviour changes relative to the pre-refactor version
------------------------------------------------------
* ppm axes are nmrglue's (Bruker's ``OFFSET - i*SW_p/(SF*SI)``) instead of
  ``linspace(OFFSET, OFFSET - SW, FTSIZE)``; integrals over a peak change
  by up to ~1 % and positions by up to ~0.5 ppm at the far edge.
* ``set_fid`` removes Bruker's digital filter with nmrglue (also works with
  a fractional ``GRPDLY``), and its time axis uses the true spacing of
  complex points, ``1/SW_h`` (it was half of that).
* ``Espectro.spec`` is always ``real + 1j*imag``, so crops can no longer
  leave it out of sync.
"""

# Names re-exported to the scripts through ``from Datos import *``.
import numpy as np
from Acqus import *
from Procs import *
from Espectro import *
from Fid import *
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.stats import linregress
from pathlib import Path
import os
import nmrglue as ng

import nmrtools as nt


def _part_name(absolute=False, imaginary=False):
    """
    Translate the legacy ``absolute``/``imaginary`` flags to an nmrtools part.

    Parameters
    ----------
    absolute : bool
        Use ``|spec|``.
    imaginary : bool
        Use the imaginary part (ignored when ``absolute`` is set).

    Returns
    -------
    str
        ``"abs"``, ``"imag"`` or ``"real"``.
    """
    if absolute:
        return "abs"
    if imaginary:
        return "imag"
    return "real"


class Datos(object):
    """
    Bruker experiment: acquisition parameters and, optionally, the raw FID.

    Parameters
    ----------
    directorio : str
        Experiment folder (the one containing ``acqus``).
    set_fid : bool, optional
        Read the raw FID into :attr:`fid` (see :meth:`set_fid`).
    read_pp : bool, optional
        Read the pulse-program description into :attr:`pulseprogram`.
    ndim : int, optional
        Kept for compatibility; the dimensionality is read from the data.

    Attributes
    ----------
    directorio : str
        Experiment folder.
    acqus : Acqus
        Acquisition parameters (``acqus.dic`` holds every parameter).
    pulseprogram : str or None
        Second line of the ``pulseprogram`` file (its description).
    title : str
        TopSpin title of ``pdata/1``.
    fid : Fid
        Raw data container, filled by :meth:`set_fid`.
    nucleo : str
        Observed nucleus (``NUC1``).
    gamma : float
        Gyromagnetic ratio of :attr:`nucleo`, in rad s^-1 T^-1.
    """

    def __init__(self, directorio, set_fid=False,
                 read_pp=True, ndim=1):

        self.directorio = directorio
        self.acqus = Acqus(directorio)
        self.pulseprogram = None
        self.ndim = ndim
        if read_pp:
            self.pulseprogram = self.getPulseProg()
        try:
            self.title = self.get_title()
        except FileNotFoundError:
            self.title = "No title file"
        self.fid = Fid()

        self.nucleo = self.acqus.dic["NUC1"]
        self.gamma = self.set_gamma()  # rad/(s*T)
        if set_fid:
            self.set_fid()

    def get_title(self, p_dir=1):
        """
        Read the TopSpin title of a processing folder.

        Parameters
        ----------
        p_dir : int, optional
            Processing number. Default 1.

        Returns
        -------
        str
            Contents of ``pdata/<p_dir>/title``.
        """
        title_file = os.path.join(self.directorio, 'pdata', str(p_dir), 'title')
        with open(title_file, 'r') as file:
            data = file.read()
        return data

    def getPulseProg(self):
        """
        Description line of the pulse program used.

        Returns
        -------
        str
            Second line of the ``pulseprogram`` file, without the leading
            ``;`` and the trailing newline.

        Raises
        ------
        FileNotFoundError
            If there is no file with ``pulseprogram`` in its name.
        """
        files = list(Path(self.directorio).glob("*pulseprogram*"))
        if not files:
            raise FileNotFoundError(f"No se encontró ningún archivo que contenga 'pulseprogram' en {self.directorio}")

        exact_match = [f for f in files if f.name == "pulseprogram"]
        if len(exact_match) == 1:
            ppfile = str(exact_match[0])
        elif len(exact_match) > 1:
            print("Warning: Se encontraron múltiples archivos llamados exactamente 'pulseprogram'. Se usará el primero.")
            ppfile = str(exact_match[0])
        else:
            if len(files) > 1:
                print("Warning: Se encontraron varios archivos que contienen 'pulseprogram', pero ninguno coincide exactamente.")
            ppfile = str(files[0])

        with open(ppfile, 'r') as file:
            data = file.readlines()[1]
            data = data[1:-1]
        return data

    def set_gamma(self):
        """
        Gyromagnetic ratio of the observed nucleus.

        Returns
        -------
        float
            γ in rad s^-1 T^-1, from :data:`nmrtools.GAMMA`.

        Raises
        ------
        KeyError
            If the nucleus is not in :data:`nmrtools.GAMMA`.
        """
        return nt.gamma(self.nucleo)

    def set_fid(self):
        """
        Read the raw FID (or ``ser``) into :attr:`fid` with :func:`nmrtools.load_fid`.

        The digital-filter group delay is removed with nmrglue, intensities
        are divided by NS·RG and the time axis has spacing ``1/SW_h``.
        Pseudo-2D data gives a 2D signal (one FID per row).
        """
        raw = nt.load_fid(self.directorio)
        self.fid.set_timeAxis(raw.time)
        self.fid.set_signal(raw.data)
        self.fid.set_size(raw.data.shape)
        self.fid.set_ppm(self.acqus.SFO1)


# ----------------------------HERENCIAS-----------------------------------------


class DatosProcesados(Datos):
    """
    Processed 1D spectrum (``pdata/<p_dir>``).

    Parameters
    ----------
    directorio : str
        Experiment folder.
    p_dir : int, optional
        Processing number. Default 1.
    ppmRange : sequence of float, optional
        Default integration range for :meth:`Integrar`.
    read_pp : bool, optional
        Read the pulse-program description.
    nomrmalize_NS : bool, optional
        Divide by NS as well as RG (default). The misspelt name is kept for
        compatibility with existing scripts.
    ndim : int, optional
        Kept for compatibility.

    Attributes
    ----------
    dataset : nmrtools.Dataset
        The underlying immutable dataset, as loaded (before any
        ``espectro.ppmSelect``).
    espectro : Espectro
        Mutable spectrum used by the legacy methods.
    procs : Procs
        Processing parameters.
    """

    def __init__(self, directorio, p_dir=1,
                 ppmRange=None,
                 read_pp=False,
                 nomrmalize_NS=True,
                 ndim=1):
        Datos.__init__(self, directorio, read_pp=read_pp, ndim=ndim)
        self.p_dir = p_dir
        try:
            self.title = self.get_title(p_dir)
        except FileNotFoundError:
            self.title = "No title file"
        self.procs = Procs(directorio, p_dir)
        self.espectro = Espectro()
        self.ppmRange = ppmRange
        self.normalize_NS = nomrmalize_NS
        self.dataset = None
        self.set_espectro()

    def set_espectro(self):
        """
        Load the processed spectrum with :func:`nmrtools.load` into :attr:`espectro`.
        """
        normalize = "ns_rg" if self.normalize_NS else "rg"
        ds = nt.load(self.directorio, procno=self.p_dir, normalize=normalize)
        self.dataset = ds
        self.espectro.set_real(ds.real.copy())
        self.espectro.set_imag(ds.imag.copy())
        self.espectro.set_size(ds.spec.shape)
        self.espectro.set_ppmAxis(ds.ppm.copy())

    def crear_ppmAxis(self, procs, acqus=None, SF=None):
        """
        ppm axis of one dimension from its processing parameters.

        Uses Bruker's definition, identical to nmrglue's axis.

        Parameters
        ----------
        procs : Procs
            Processing parameters of the dimension (``procs`` or ``proc2s``).
        acqus : Acqus, optional
            Unused; kept for compatibility.
        SF : float, optional
            Unused; kept for compatibility.

        Returns
        -------
        numpy.ndarray
            ``OFFSET - i * SW_p / (SF * SI)`` for ``i = 0 .. SI-1``.
        """
        d = procs.dic
        return d["OFFSET"] - np.arange(d["SI"]) * d["SW_p"] / d["SF"] / d["SI"]

    def _range_or_default(self, ppmRange):
        """
        Resolve the integration range: argument, then :attr:`ppmRange`, then all.

        Parameters
        ----------
        ppmRange : sequence of float or None
            Range passed by the caller.

        Returns
        -------
        sequence of float or None
            ``None`` means the whole axis.
        """
        if ppmRange is None:
            if self.ppmRange is not None:
                return self.ppmRange
            print("WARNING: No se ha definido rango de integracion")
        return ppmRange

    def Integrar(self, ppmRange=None, absolute=False, imaginary=False):
        """
        Integral of the (possibly cropped) spectrum over a ppm range.

        Parameters
        ----------
        ppmRange : sequence of float, optional
            Limits; defaults to :attr:`ppmRange`, then the whole axis.
        absolute : bool, optional
            Integrate ``|spec|``.
        imaginary : bool, optional
            Integrate the imaginary part.

        Returns
        -------
        numpy.ndarray
            The integral (0-d for 1D, one value per row for 2D). Also
            stored in :attr:`signal`.
        """
        ppmRange = self._range_or_default(ppmRange)
        signal = np.array(nt.integrate(self.espectro.ppmAxis, self.espectro.spec,
                                       ppmRange, _part_name(absolute, imaginary)))
        self.signal = signal
        return signal

    def Mean_ppm(self, ppmRange=None, absolute=False):
        r"""
        Intensity-weighted mean chemical shift of a region.

        ``mean = \int ppm S(ppm) dppm / \int S(ppm) dppm``

        Parameters
        ----------
        ppmRange : sequence of float, optional
            Limits; defaults to :attr:`ppmRange`, then the whole axis.
        absolute : bool, optional
            Use ``|spec|`` as weight.

        Returns
        -------
        numpy.ndarray
            Mean position in ppm.
        """
        ppmRange = self._range_or_default(ppmRange)
        return np.array(nt.mean_ppm(self.espectro.ppmAxis, self.espectro.spec,
                                    ppmRange, _part_name(absolute)))


class DatosProcesados2D(DatosProcesados):
    """
    Processed 2D or pseudo-2D spectrum (T1, diffusion, EXSY, in-situ series...).

    Parameters
    ----------
    directorio : str
        Experiment folder.
    p_dir : int, optional
        Processing number. Default 1.
    ppmRange : sequence of float, optional
        Default integration range for :meth:`Integrar2D`.
    read_pp : bool, optional
        Read the pulse-program description.
    extract_vdfactor : bool, optional
        Look for ``vd*N`` in the pulse program when building delays.
    ndim : int, optional
        Kept for compatibility.
    remove_empty : bool, optional
        Drop trailing rows that were never acquired.

    Attributes
    ----------
    signal : numpy.ndarray or None
        One integral per row, after :meth:`Integrar2D`.
    proc2s, acqu2s : Procs, Acqus
        Parameters of the indirect dimension.
    """

    def __init__(self, directorio, p_dir=1, ppmRange=None,
                 read_pp=False, extract_vdfactor=False, ndim=2, remove_empty=False):
        DatosProcesados.__init__(self, directorio, p_dir=p_dir,
                                 read_pp=read_pp, ndim=ndim)
        self.ndim = 2
        self.proc2s = Procs(directorio, p_dir, dim2=True)
        self.acqu2s = Acqus(directorio, dim2=True)
        self.espectro.set_ppmAxisInd(self.dataset.ppm_ind.copy())
        self.espectro.crear_ppmGrid()
        self.signal = None
        self.extract_vdfactor = extract_vdfactor

        self.ppmRange = ppmRange
        self.remove_empty = remove_empty
        if remove_empty:
            self._remove_empty_rows()

    def _remove_empty_rows(self, threshold=0.0):
        """
        Drop trailing rows whose ``|spec|`` never exceeds ``threshold``.

        Parameters
        ----------
        threshold : float, optional
            Rows with every value ``<= threshold`` count as empty.
        """
        filled = np.any(np.abs(self.espectro.spec) > threshold, axis=1)
        if not np.any(filled):
            print("WARNING: all rows are empty!")
            return
        last_idx = int(np.where(filled)[0][-1] + 1)
        self.espectro.real = self.espectro.real[:last_idx, :]
        self.espectro.imag = self.espectro.imag[:last_idx, :]
        self.espectro.ppmAxisInd = self.espectro.ppmAxisInd[:last_idx]
        self.proc2s.SI = last_idx
        self.espectro.crear_ppmGrid()

    def Integrar2D(self, ppmRange=None, absolute=False, imaginary=False):
        """
        Integral of every row over a ppm range.

        Parameters
        ----------
        ppmRange : sequence of float, optional
            Limits; defaults to :attr:`ppmRange`, then the whole axis.
        absolute : bool, optional
            Integrate ``|spec|``.
        imaginary : bool, optional
            Integrate the imaginary part.

        Returns
        -------
        numpy.ndarray
            One integral per row. Also stored in :attr:`signal`.
        """
        return self.Integrar(ppmRange, absolute=absolute, imaginary=imaginary)

    def get_vdlist(self, units='ms'):
        """
        Read the variable-delay list with :func:`nmrtools.read_vdlist`.

        Parameters
        ----------
        units : {"ms", "s"}, optional
            Units of the returned delays. Default ``"ms"``.

        Returns
        -------
        numpy.ndarray
            Delays from ``vdlist`` (or ``vdlist_path`` if set).
        """
        vdlist_path = getattr(self, "vdlist_path", None)
        vdlist = nt.read_vdlist(self.directorio, vdlist_path or "vdlist")
        if vdlist_path is not None:
            print(f"vdlist leido de {self.directorio}/{vdlist_path}")
        if units == 'ms':
            vdlist = vdlist * 1000
        return vdlist


class DatosProcesadosT1(DatosProcesados2D):
    """
    Saturation/inversion-recovery pseudo-2D with a ``vdlist``.

    Parameters
    ----------
    directorio : str
        Experiment folder.
    p_dir : int, optional
        Processing number. Default 1.
    vdlist_path : str, optional
        Delay-list file relative to ``directorio`` (default ``vdlist``).
    ndim : int, optional
        Kept for compatibility.

    Attributes
    ----------
    tau : numpy.ndarray
        Delays in ms (vdlist times the vd factor).
    signal : numpy.ndarray
        Integral of each row (whole axis until :meth:`get_T1data` is called
        with a range).
    fit : nmrtools.FitResult or None
        Last fit made by :meth:`T1fit` or :meth:`T2fit`.
    """

    def __init__(self, directorio, p_dir=1, vdlist_path=None, ndim=2):
        DatosProcesados2D.__init__(self, directorio, p_dir=p_dir, ndim=ndim)

        self.tau = None
        self.signal = None
        self.factor_vd = None
        self.vdlist_path = vdlist_path

        self.Crear_tau(directorio, p_dir)
        self.Integrar2D()

        self.tau_fit = None
        self.signal_fit = None
        self.T1params = None
        self.T1params_err = None
        self.fit = None

    def Crear_tau(self, directorio, p_dir):
        """
        Build :attr:`tau` (ms) from the vdlist and, optionally, the vd factor.

        Parameters
        ----------
        directorio : str
            Experiment folder.
        p_dir : int
            Processing number (unused; kept for compatibility).

        Returns
        -------
        int
            0, for compatibility.
        """
        vdlist = self.get_vdlist()
        factor_vd = 1
        if self.extract_vdfactor:
            factor_vd = nt.vd_factor(directorio, default=None)
            if factor_vd is None:
                print("WARNING: No se ha encontrado el factor vd. Seteando factor_vd = 1")
                factor_vd = 1
        self.factor_vd = factor_vd
        self.tau = vdlist * factor_vd
        return 0

    def get_T1data(self, ppmRange=None):
        """
        Delays and row integrals, optionally over a new ppm range.

        Parameters
        ----------
        ppmRange : sequence of float, optional
            New integration range (overwrites :attr:`ppmRange`).

        Returns
        -------
        tau : numpy.ndarray
            Delays in ms.
        signal : numpy.ndarray
            Integrals, truncated to the same length as ``tau``.
        """
        if ppmRange is not None:
            self.ppmRange = ppmRange
            self.Integrar2D()
        size = min(self.tau.size, self.signal.size)
        self.tau = self.tau[:size]
        self.signal = self.signal[:size]
        return self.tau, self.signal

    def _store_fit(self, fit, tau_fit, labels, title):
        """
        Save a fit in the legacy attributes and print it.

        Parameters
        ----------
        fit : nmrtools.FitResult
            The fit.
        tau_fit : numpy.ndarray
            Points where :attr:`signal_fit` is evaluated.
        labels : list of str
            Legacy parameter names, in fit order, for printing.
        title : str
            Header of the printed summary.

        Returns
        -------
        numpy.ndarray
            Residuals in the original order of :attr:`tau`.
        """
        self.fit = fit
        self.tau_fit = tau_fit
        self.signal_fit = fit.predict(tau_fit)
        self.T1params = np.array(list(fit.params.values()))
        self.T1stderr = np.array(list(fit.stderr.values()))
        self.T1params_err = 1.96 * self.T1stderr
        self.R2 = fit.r2
        print(title)
        for name, val, err in zip(labels, self.T1params, self.T1params_err):
            print(f"  {name} = {val:.2f} ± {err:.2f}")
        print(f"  R^2 = {fit.r2:.6f}")
        return self.signal - fit.predict(self.tau)

    def T1fit(self, model='mono'):
        """
        Exponential recovery fit with :func:`nmrtools.exp_fit`.

        Parameters
        ----------
        model : {"mono", "bi"}, optional
            ``y0 + A exp(-tau/T1)`` or the sum of two exponentials.

        Returns
        -------
        tau_fit : numpy.ndarray
            Log-spaced delays for plotting the fit.
        signal_fit : numpy.ndarray
            Fitted curve at ``tau_fit``.
        residuals : numpy.ndarray
            ``signal - fit`` at each delay.
        """
        if 'mono' in model.lower():
            kind, labels = "mono", ['A', 'T1', 'y0']
        elif 'bi' in model.lower():
            kind, labels = "bi", ['A1', 'T11', 'A2', 'T12', 'y0']
        else:
            raise ValueError("model debe ser 'mono' o 'bi'")
        fit = nt.exp_fit(self.tau, self.signal, kind, "recovery")
        tau = np.sort(self.tau)
        if tau[0] > 0:
            tau_fit = np.logspace(np.log10(tau[0]), np.log10(tau[-1]), 10 * tau.size)
        else:
            tau_fit = np.linspace(tau[0], tau[-1], 10 * tau.size)
        residuals = self._store_fit(fit, tau_fit, labels, f"Ajuste {model}-exponencial de T1:")
        return self.tau_fit, self.signal_fit, residuals

    def T2fit(self):
        """
        Mono-exponential decay fit with :func:`nmrtools.exp_fit`.

        Returns
        -------
        tau_fit : numpy.ndarray
            512 delays from 0 to the last delay.
        signal_fit : numpy.ndarray
            Fitted curve at ``tau_fit``.
        residuals : numpy.ndarray
            ``signal - fit`` at each delay.
        """
        fit = nt.exp_fit(self.tau, self.signal, "mono", "decay")
        tau_fit = np.linspace(0, self.tau[-1], 512)
        residuals = self._store_fit(fit, tau_fit, ['A', 'T2', 'y0'], "Ajuste exponencial de T2:")
        return self.tau_fit, self.signal_fit, residuals


# ------------------------------------------------------------------------------
class DatosProcesadosDiff(DatosProcesados2D):
    """
    Bipolar STE diffusion pseudo-2D.

    Parameters
    ----------
    directorio : str
        Experiment folder.
    p_dir : int, optional
        Processing number. Default 1.
    factor_b : float, optional
        Gradient calibration factor applied to the b-values.
    bmax, bmin : float, optional
        b-range kept by :meth:`Recortar_datos` (units of 1e9 s/m^2).
    gpshape : {"sin", "rect"}, optional
        Gradient pulse shape.

    Attributes
    ----------
    bvalue : numpy.ndarray
        b-values in units of 1e9 s/m^2 (first point dropped).
    signal : numpy.ndarray
        Row integrals (S vs b).
    fit_results : list
        ``[D, uD, r_squared]`` with D in 1e-9 m^2/s, after :meth:`Diff1_fit`.
    """

    def __init__(self, directorio, p_dir=1, factor_b=1, bmax=np.inf, bmin=0, gpshape='sin'):
        DatosProcesados2D.__init__(self, directorio, p_dir=p_dir)

        self.bvalue = None
        self.signal = None
        self.factor_b = factor_b
        self.bmax = bmax
        self.bmin = bmin
        self.absolute = None
        # For bipolar STE, delta is 2*P30 and bigDelta is D20.
        self.delta = 2*self.acqus.P[30]*1e-3  # ms (P30 is in us)
        self.bigDelta = self.acqus.D[20]*1e3  # ms (D20 is in s)
        self.gpmax = self.acqus.GPZ[6]
        self.gpvalue = None

        self.Crear_bvalue(gpshape)
        # drop the first point:
        Npts = self.acqu2s.TD
        self.espectro.IndirectDimensionSelect([1, Npts])
        self.bvalue = self.bvalue[1:]

        self.Integrar2D()

        self.signal_fit = None

    def Crear_bvalue(self, gpshape):
        """
        Build :attr:`bvalue` with :func:`nmrtools.bvalues`.

        Parameters
        ----------
        gpshape : str
            Contains ``"sin"`` or ``"rect"`` (case-insensitive).

        Returns
        -------
        int
            0, for compatibility.

        Raises
        ------
        ValueError
            If ``gpshape`` names neither shape.
        """
        if 'sin' in gpshape.lower():
            shape = 'sin'
        elif 'rect' in gpshape.lower():
            shape = 'rect'
        else:
            raise ValueError("gpshape debe ser 'sin' o 'rect'")
        gplist = np.linspace(0, self.acqus.gp, self.acqu2s.TD)
        self.gpvalue = gplist[1:]
        b = nt.bvalues(gplist, self.delta*1e-3, self.bigDelta*1e-3, self.gamma,
                       g_max=12, shape=shape)
        self.bvalue = b * 1e-9 * self.factor_b
        return 0

    def get_Diffdata(self, ppmRange=None, absolute=False):
        """
        b-values and integrals, optionally over a new ppm range.

        Parameters
        ----------
        ppmRange : sequence of float, optional
            New integration range (overwrites :attr:`ppmRange`).
        absolute : bool, optional
            Integrate ``|spec|``.

        Returns
        -------
        bvalue : numpy.ndarray
        signal : numpy.ndarray
        """
        self.absolute = absolute
        if ppmRange is not None:
            self.ppmRange = ppmRange
            self.Integrar2D(absolute=absolute)
        self.Recortar_datos()
        return self.bvalue, self.signal

    def Diff1_fit(self):
        """
        Linear fit of ``ln S`` vs b with :func:`nmrtools.diffusion_fit`.

        Returns
        -------
        bvalue : numpy.ndarray
        signal_fit : numpy.ndarray
            Fitted attenuation at each b.
        residuals : numpy.ndarray
            ``signal - signal_fit``.
        """
        fit = nt.diffusion_fit(self.bvalue, self.signal)
        signal_fit = fit.predict(self.bvalue)
        residuals = self.signal - signal_fit
        self.fit_results = [fit.params["D"], fit.stderr["D"], fit.r2]
        print(f"Ajuste lineal de Difusion:\n"
              f"  D =  {fit.params['D']:.8f} 10^-9 m^2/s\n"
              f"  Rsquared = {fit.r2:.6f}")
        self.signal_fit = signal_fit
        return self.bvalue, self.signal_fit, residuals

    def Recortar_datos(self):
        """
        Keep only the points with ``bmin <= b <= bmax``.

        Returns
        -------
        int
            0, for compatibility.
        """
        if self.bmax == np.inf:
            bmax = self.bvalue[-1]
        else:
            bmax = self.bmax
        bmin = self.bmin
        if bmax != np.inf or bmin != 0:
            bvalue = self.bvalue
            signal = self.signal

            bmed = (bmax+bmin)/2.00
            semiancho = (bmax-bmin)/2.00
            condicion = np.abs(bvalue-bmed) <= semiancho
            self.bvalue = bvalue[condicion]
            self.gpvalue = self.gpvalue[condicion]
            self.signal = signal[condicion]
        return 0
