# -*- coding: utf-8 -*-
"""
Relaxation analysis (T1 recovery or T2 decay) of one or many pseudo-2D experiments.

Replaces the family of ``read_T1data*.py`` / ``read_pseudo2Ddata_integrate_
DNP_LiMetalT1/T2*.py`` scripts with one script built on :mod:`nmrtools`.
For every experiment in :data:`EXPERIMENTS` it

1. loads the processed pseudo-2D and its delay list (any vdlist file),
2. integrates the spectra over fixed ppm ranges, or over a window centred on
   the spectral maximum,
3. fits each integral against the delays with a mono- or bi-exponential
   recovery (T1) or decay (T2),
4. draws spectra, fits and residuals.

Results from all experiments are gathered in one table, ``summary.csv``
(one row per experiment and range), plotted against any column you choose,
e.g. temperature or start time. That table is the input for follow-up
analysis such as temperature calibrations or Arrhenius fits.

Files written to :data:`SAVEPATH` when :data:`SAVE` is true
-----------------------------------------------------------
``summary.csv``
    One row per experiment and range: label, expno, start time, your extra
    experiment keys, ppm range, ppm of maximum, and the fit (``T1_s``,
    ``T1_err_s`` = 1 standard error, ``T1_ci95_s``, ``A``, ``y0``, ``r2``;
    for bi-exponential fits also short/long components and fractions).
    Times are in seconds. For ``KIND = "decay"`` the columns say T2.
``summary.png``
    The relaxation time of every range against :data:`SUMMARY_X`.
``<label>_T1.dat``
    Delays (s) and the integral of every range, tab-separated.
``<label>_lastSpectrum.dat``
    The reference spectrum (ppm, real, imag) used to place centred windows.
``<label>_T1.png``
    Spectra coloured by delay, integrals with fits, residuals.
``run_record.json``
    The parameters of this run and the git commit (see :mod:`provenance`).

Run it as a whole, or cell by cell (``#%%``) in Spyder/VS Code; the results
stay available in :data:`results` and :data:`summary`.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import nmrtools as nt
from paths import NMRDATA, HOME
from provenance import write_run_record

#%% ========================== PARAMETERS ==================================
#: Folder that contains the experiment numbers.
PATH = rf"{NMRDATA}\500\2025-08-05_PEO-PTT-solid-electrolyte_VT"

#: Output folder for all experiments of this run.
SAVEPATH = rf"{HOME}\PolymerElectrolyte\Analysis\2025-08_500MHz_VT\T1_79Br"

#: ``"recovery"`` for T1 (saturation or inversion recovery) or ``"decay"``
#: for T2 / CPMG-like decays.
KIND = "recovery"

#: ``"mono"`` or ``"bi"`` exponential.
MODEL = "mono"

#: Delay list, relative to each experiment folder (e.g. "lists/vd/sm2974").
VDLIST = "vdlist"

#: Multiplier of the listed delays, or ``"auto"`` to read ``vd*N`` from the
#: pulse program.
VD_FACTOR = 1

#: ppm region shown in the plots and searched for the spectral maximum.
PLOT_RANGE = (100, -100)

#: Signal definition, one of:
#:
#: * ``("centred", width)``: one window of ``width`` ppm centred on the
#:   maximum of the reference spectrum (last delay for recovery, first for
#:   decay) inside PLOT_RANGE;
#: * ``("fixed", [(high, low), ...])``: one or more fixed ppm windows.
SIGNAL = ("centred", 20)

#: One dict per experiment. ``label`` and ``expno`` are required. Optional
#: keys ``vdlist``, ``vd_factor`` and ``signal`` override the defaults above
#: for that experiment. Any other key (e.g. a temperature) is copied to
#: ``summary.csv`` and can be used as SUMMARY_X.
T_NOMINAL = [-20, -10, 0, 10, 20, 30, 40]
EXPERIMENTS = [dict(label=f"T{t:+d}C", expno=20 + 10 * i, T_nominal_C=t)
               for i, t in enumerate(T_NOMINAL)]

#: Column of summary.csv used as x axis of summary.png (e.g. "T_nominal_C",
#: "start_time", "expno" or "label").
SUMMARY_X = "T_nominal_C"

#: Write tables, figures and the run record.
SAVE = True


#%% ========================== ANALYSIS ====================================
QUANTITY = "T1" if KIND == "recovery" else "T2"


def delays(expdir, vdlist=VDLIST, vd_factor=VD_FACTOR):
    """
    Delays of one experiment, in seconds.

    Parameters
    ----------
    expdir : str
        Experiment folder.
    vdlist : str, optional
        Delay-list file relative to ``expdir``.
    vd_factor : float or "auto", optional
        Multiplier, or ``"auto"`` to read it from the pulse program.

    Returns
    -------
    numpy.ndarray
    """
    factor = nt.vd_factor(expdir) if vd_factor == "auto" else float(vd_factor)
    return nt.read_vdlist(expdir, vdlist) * factor


def signal_windows(ppm, spec, signal, kind=KIND, search=PLOT_RANGE):
    """
    ppm windows to integrate, and the position of the spectral maximum.

    Parameters
    ----------
    ppm : numpy.ndarray
        ppm axis.
    spec : numpy.ndarray
        Spectra, one row per delay.
    signal : tuple
        ``("centred", width)`` or ``("fixed", [(high, low), ...])``.
    kind : {"recovery", "decay"}, optional
        Picks the reference spectrum: last row for recovery, first for decay.
    search : tuple of float, optional
        ppm region where the maximum is searched.

    Returns
    -------
    windows : list of tuple
        ``(high, low)`` ppm pairs.
    ppm_of_max : float
        Position of the maximum of the reference spectrum within ``search``.
    reference : numpy.ndarray
        The reference spectrum (complex).
    """
    reference = spec[-1] if kind == "recovery" else spec[0]
    lo, hi = min(search), max(search)
    inside = (ppm >= lo) & (ppm <= hi)
    ppm_of_max = float(ppm[inside][np.argmax(reference.real[inside])])
    mode, value = signal
    if mode == "centred":
        windows = [(ppm_of_max + value / 2, ppm_of_max - value / 2)]
    elif mode == "fixed":
        windows = [tuple(w) for w in value]
    else:
        raise ValueError("signal must be ('centred', width) or ('fixed', [ranges])")
    return windows, ppm_of_max, reference


def fit_row(fit, model):
    """
    Flatten a fit into summary columns.

    Parameters
    ----------
    fit : nmrtools.FitResult or None
        The fit; ``None`` if it failed (gives NaN columns).
    model : {"mono", "bi"}
        Model that was fitted.

    Returns
    -------
    dict
        ``<Q>_s``, ``<Q>_err_s``, ``<Q>_ci95_s``, ``A``, ``y0``, ``r2`` for
        mono; for bi the short/long times, errors, amplitudes and
        fractions, ``y0`` and ``r2``. ``<Q>`` is T1 or T2.
    """
    q = QUANTITY
    if model == "mono":
        if fit is None:
            return {f"{q}_s": np.nan, f"{q}_err_s": np.nan, f"{q}_ci95_s": np.nan,
                    "A": np.nan, "y0": np.nan, "r2": np.nan}
        return {f"{q}_s": fit.params["tau"], f"{q}_err_s": fit.stderr["tau"],
                f"{q}_ci95_s": fit.ci95["tau"], "A": fit.params["A"],
                "y0": fit.params["y0"], "r2": fit.r2}
    keys = [f"{q}_short_s", f"{q}_short_err_s", "A_short", "frac_short",
            f"{q}_long_s", f"{q}_long_err_s", "A_long", "frac_long", "y0", "r2"]
    if fit is None:
        return dict.fromkeys(keys, np.nan)
    p, e = fit.params, fit.stderr
    comps = sorted([(p["tau1"], e["tau1"], p["A1"]), (p["tau2"], e["tau2"], p["A2"])])
    total = comps[0][2] + comps[1][2]
    return {f"{q}_short_s": comps[0][0], f"{q}_short_err_s": comps[0][1],
            "A_short": comps[0][2], "frac_short": comps[0][2] / total,
            f"{q}_long_s": comps[1][0], f"{q}_long_err_s": comps[1][1],
            "A_long": comps[1][2], "frac_long": comps[1][2] / total,
            "y0": p["y0"], "r2": fit.r2}


def analyse_experiment(exp):
    """
    Load, integrate and fit one experiment.

    Parameters
    ----------
    exp : dict
        One entry of :data:`EXPERIMENTS`.

    Returns
    -------
    dict
        ``ds`` (the dataset), ``tau`` (s), ``windows``, ``ppm_of_max``,
        ``reference`` spectrum, ``signals`` (one array per window),
        ``fits`` (one :class:`nmrtools.FitResult` or ``None`` per window).
    """
    expdir = os.path.join(PATH, str(exp["expno"]))
    ds = nt.load(expdir)
    tau = delays(expdir, exp.get("vdlist", VDLIST), exp.get("vd_factor", VD_FACTOR))
    n = min(tau.size, ds.spec.shape[0])   # TopSpin may pad the 2D with empty rows
    tau, spec = tau[:n], ds.spec[:n]

    windows, ppm_of_max, reference = signal_windows(ds.ppm, spec, exp.get("signal", SIGNAL))
    signals, fits = [], []
    for window in windows:
        y = nt.integrate(ds.ppm, spec, window)
        signals.append(y)
        try:
            fits.append(nt.exp_fit(tau, y, MODEL, KIND))
        except RuntimeError as err:
            print(f"[{exp['label']}] fit failed for {window} ppm: {err}")
            fits.append(None)
    return dict(ds=ds, spec=spec, tau=tau, windows=windows, ppm_of_max=ppm_of_max,
                reference=reference, signals=signals, fits=fits)


def summary_rows(exp, result):
    """
    Rows of ``summary.csv`` for one experiment (one per window).

    Parameters
    ----------
    exp : dict
        The experiment entry (extra keys are copied to the row).
    result : dict
        Output of :func:`analyse_experiment`.

    Returns
    -------
    list of dict
    """
    extra = {k: v for k, v in exp.items() if k not in ("label", "expno", "vdlist", "vd_factor", "signal")}
    rows = []
    for i, (window, fit) in enumerate(zip(result["windows"], result["fits"])):
        rows.append({"label": exp["label"], "expno": exp["expno"],
                     "start_time": result["ds"].start_time, **extra,
                     "window": i, "ppm_high": max(window), "ppm_low": min(window),
                     "ppm_of_max": result["ppm_of_max"], "model": MODEL,
                     **fit_row(fit, MODEL)})
    return rows


#%% ========================== FIGURES =====================================
def plot_experiment(result, title):
    """
    Spectra coloured by delay, integrals with fits, and residuals.

    Parameters
    ----------
    result : dict
        Output of :func:`analyse_experiment`.
    title : str
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    ds, tau = result["ds"], result["tau"]
    fig, (ax_spec, ax_fit, ax_res) = plt.subplots(1, 3, figsize=(14, 4))
    colors = plt.get_cmap("viridis")(np.linspace(0, 1, tau.size))
    ppm, spec = nt.crop(ds.ppm, result["spec"], PLOT_RANGE)
    for k in range(tau.size):
        ax_spec.plot(ppm, spec[k].real, color=colors[k], lw=0.8)
    for window in result["windows"]:
        ax_spec.axvspan(min(window), max(window), color="grey", alpha=0.2)
    ax_spec.set(xlim=(max(PLOT_RANGE), min(PLOT_RANGE)), xlabel=r"$\delta$ (ppm)",
                ylabel="NMR signal [a.u.]", title="spectra (colour: delay)")

    for i, (y, fit) in enumerate(zip(result["signals"], result["fits"])):
        window = result["windows"][i]
        line, = ax_fit.plot(tau, y, "o", label=f"{max(window):.4g} to {min(window):.4g} ppm")
        if fit is None:
            continue
        x, curve = fit.curve()
        ax_fit.plot(x, curve, "-", color=line.get_color())
        ax_res.plot(fit.x, fit.residuals, "o-", color=line.get_color())
        row = fit_row(fit, MODEL)
        if MODEL == "mono":
            text = f"{QUANTITY} = {row[f'{QUANTITY}_s']:.4g} ± {row[f'{QUANTITY}_ci95_s']:.2g} s"
        else:
            text = (f"{QUANTITY} = {row[f'{QUANTITY}_short_s']:.3g} s ({row['frac_short']:.0%}), "
                    f"{row[f'{QUANTITY}_long_s']:.3g} s ({row['frac_long']:.0%})")
        ax_fit.text(0.02, 0.97 - 0.07 * i, text, transform=ax_fit.transAxes,
                    va="top", fontsize=8, color=line.get_color())
    for ax in (ax_fit, ax_res):
        ax.set_xscale("log")
        ax.set_xlabel("delay (s)")
        ax.grid(alpha=.3)
    ax_fit.set(ylabel="integral [a.u.]", title=f"{MODEL}-exponential {KIND}")
    ax_fit.legend(fontsize=7, loc="lower right" if KIND == "recovery" else "upper right")
    ax_res.axhline(0, color="k", lw=0.5)
    ax_res.set(ylabel="residual", title="residuals")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_summary(table, x=SUMMARY_X):
    """
    Relaxation time of every window against a column of the summary table.

    Parameters
    ----------
    table : pandas.DataFrame
        The summary table.
    x : str, optional
        Column used as x axis.

    Returns
    -------
    matplotlib.figure.Figure
    """
    q = QUANTITY
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for window, part in table.groupby("window"):
        label = f"{part.ppm_high.iloc[0]:.4g} to {part.ppm_low.iloc[0]:.4g} ppm"
        if MODEL == "mono":
            ax.errorbar(part[x], part[f"{q}_s"], yerr=part[f"{q}_ci95_s"], fmt="o-",
                        capsize=3, label=label)
        else:
            ax.plot(part[x], part[f"{q}_short_s"], "o-", label=f"{label}, short")
            ax.plot(part[x], part[f"{q}_long_s"], "s--", label=f"{label}, long")
    ax.set(xlabel=x, ylabel=f"{q} (s)", yscale="log")
    ax.grid(alpha=.3, which="both")
    if table.window.nunique() > 1 or MODEL == "bi":
        ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


#%% ========================== OUTPUT ======================================
def save_experiment(folder, label, result, fig):
    """
    Write the integrals, the reference spectrum and the figure of one experiment.

    Parameters
    ----------
    folder : str
        Output folder.
    label : str
        Experiment label, used as file prefix.
    result : dict
        Output of :func:`analyse_experiment`.
    fig : matplotlib.figure.Figure
        Figure from :func:`plot_experiment`.
    """
    header = "tau [s]\t" + "\t".join(f"[{max(w):.6g}, {min(w):.6g}] ppm" for w in result["windows"])
    np.savetxt(os.path.join(folder, f"{label}_{QUANTITY}.dat"),
               np.column_stack([result["tau"], *result["signals"]]), header=header, delimiter="\t")
    ref = result["reference"]
    np.savetxt(os.path.join(folder, f"{label}_lastSpectrum.dat"),
               np.column_stack([result["ds"].ppm, ref.real, ref.imag]),
               header="ppm\treal\timag", delimiter="\t")
    fig.savefig(os.path.join(folder, f"{label}_{QUANTITY}.png"), dpi=150)


#%% ========================== RUN =========================================
if __name__ == "__main__":
    #: Results per experiment label (see :func:`analyse_experiment`).
    results = {}
    rows = []
    if SAVE:
        os.makedirs(SAVEPATH, exist_ok=True)
    for exp in EXPERIMENTS:
        result = analyse_experiment(exp)
        results[exp["label"]] = result
        rows += summary_rows(exp, result)
        fig = plot_experiment(result, f"{exp['label']} (expno {exp['expno']})")
        if SAVE:
            save_experiment(SAVEPATH, exp["label"], result, fig)
            plt.close(fig)

    #: All fits, one row per experiment and window.
    summary = pd.DataFrame(rows)
    print(summary.drop(columns=["start_time"]).to_string(index=False))
    fig_summary = plot_summary(summary) if len(EXPERIMENTS) > 1 else None

    if SAVE:
        summary.to_csv(os.path.join(SAVEPATH, "summary.csv"), index=False)
        if fig_summary is not None:
            fig_summary.savefig(os.path.join(SAVEPATH, "summary.png"), dpi=150)
        write_run_record(SAVEPATH, dict(
            PATH=PATH, KIND=KIND, MODEL=MODEL, VDLIST=VDLIST, VD_FACTOR=VD_FACTOR,
            PLOT_RANGE=PLOT_RANGE, SIGNAL=SIGNAL, EXPERIMENTS=EXPERIMENTS,
            SUMMARY_X=SUMMARY_X), __file__)
    plt.show()
