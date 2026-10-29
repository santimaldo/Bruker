# -*- coding: utf-8 -*-
"""
Kpopt field-sweep analysis: NMR signal as a function of the bsms field offset.

Pilot rewrite of ``read_paeudo2Ddata_Kpopt_Gabriel.py`` on top of
:mod:`nmrtools`. It reads a pseudo-2D Kpopt experiment (one spectrum per
bsms value, listed in ``Klog_opt``), puts every spectrum back on a common
chemical-shift scale, and follows the integral, centre and width of a few
ppm regions as the field changes.

For each experiment in :data:`EXPERIMENTS` it writes, to
``SAVEPATH_BASE/<label>/``:

``integrals.csv``
    One row per bsms value. The first columns are those of the original
    script (``integrals``, ``bsms``, ``center``, ``width``,
    ``bsms_centred_area``, ``bsms_centred_area_mT``, ``bsms_centred``,
    ``bsms_centred_mT``) for the first ppm range, all real-valued. They are
    followed by ``area_re_i``, ``area_im_i``, ``center_i``, ``width_i`` for
    every range ``i`` of :data:`PPM_RANGES`, and ``height`` (spectrum maximum).
``integrals_results.csv``
    bsms of the maximum integral, of the minimum centre and of the maximum
    height, and the mT-per-bsms conversion factor.
``spectra/specNNN.dat``
    Each field-corrected spectrum: ppm, real, imaginary.
``2dMap_figure.png``, ``2dMap_figure_zoom.png``, ``slices.png``
    Integrals and spectral map against ΔB, and spectra at selected fields.

Run it as a whole, or cell by cell (``#%%``) in Spyder/VS Code; the results
stay available in :data:`results` for interactive use.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import nmrtools as nt
from kpopt_utils import read_bsms_field
from paths import NMRDATA_LOCAL, PROJECTS

#%% ========================== PARAMETERS ==================================
#: Folder that contains the experiment numbers.
PATH = rf"{NMRDATA_LOCAL}\400dnp\2026-09-10_3.2mm_LiMetal_LiF_RC1"

#: Results go to ``SAVEPATH_BASE/<label>/``.
SAVEPATH_BASE = rf"{PROJECTS}\LiMetal\DNP\analysis\2026-09_Limetal_LiF\0kHz"

#: One entry per Kpopt experiment:
#:
#: ``label``
#:     Sample label, also the output sub-folder.
#: ``expno``
#:     Experiment number.
#: ``ppm_offset``
#:     Subtracted from the ppm axis (chemical-shift calibration).
#: ``bsms_zero_field``
#:     bsms value taken as ΔB = 0. ``None`` uses the bsms of the tallest
#:     spectrum.
EXPERIMENTS = [
    dict(label="kpopt_p1-4.35us_p90", expno=265, ppm_offset=0.0, bsms_zero_field=-4750),
    # dict(label="kpopt_p1-1us", expno=262, ppm_offset=-98.0, bsms_zero_field=-4750),
]

#: ppm regions to follow. The first one gives the main ``integrals`` columns.
PPM_RANGES = [(300, 10), (0, -50), (300, 230), (230, 10)]

#: Regions (indices into PPM_RANGES) and markers drawn in the integral panel.
SUMMARY_RANGES = [(0, "o"), (2, "^"), (3, "x")]

#: ppm limits of the spectral plots.
PLOT_RANGE = (300, -50)

#: ΔB values (mT) whose spectra are drawn in ``slices.png``.
FIELDS_TO_PLOT_MT = [0.027, 0, -0.027, -0.083, -0.11, -0.137, -0.25, 0.15]

#: ΔB limits (mT) of ``2dMap_figure_zoom.png``.
ZOOM_MT = (0.2, -0.4)

#: Vertical guide lines (ppm) in ``slices.png``.
SLICE_GUIDES_PPM = (230, 10)

#: Colour-scale maximum of the map, as a fraction of the largest intensity.
MAP_VMAX_FRACTION = 1.0

#: Write the CSV, .dat and PNG files.
SAVE = True

# Field calibration: one bsms unit moves the lines by PPM_PER_BSMS ppm, at
# B0_TESLA. ΔB in mT follows from the relative shift.
PPM_PER_BSMS = -1.1665e-2
B0_TESLA = 9.42
MT_PER_BSMS = PPM_PER_BSMS * 1e-6 * B0_TESLA * 1e3


#%% ========================== ANALYSIS ====================================
def analyse_kpopt(ds, bsms, ppm_ranges, ppm_offset=0.0):
    """
    Field-correct a Kpopt series and measure each ppm region.

    Parameters
    ----------
    ds : nmrtools.Dataset
        Pseudo-2D Kpopt dataset, one row per bsms value.
    bsms : numpy.ndarray
        bsms value of each row (from ``Klog_opt``).
    ppm_ranges : list of tuple
        Regions to integrate, as ``(high, low)`` ppm pairs.
    ppm_offset : float, optional
        Subtracted from the ppm axis before anything else.

    Returns
    -------
    ppm : numpy.ndarray
        Calibrated ppm axis, shared by all rows.
    spec : numpy.ndarray
        Field-corrected complex spectra, shape ``(len(bsms), len(ppm))``.
        Row ``k`` is shifted by ``-PPM_PER_BSMS * bsms[k]`` ppm, in whole
        points, so that all rows share the same chemical-shift scale.
    table : pandas.DataFrame
        Per-row measurements: ``bsms``, ``height`` and, for every range
        ``i``, ``area_re_i``, ``area_im_i``, ``center_i``, ``width_i``.
    """
    ppm = ds.ppm - ppm_offset
    spec = nt.shift_ppm(ppm, ds.spec[: bsms.size], -PPM_PER_BSMS * bsms)
    table = pd.DataFrame({"bsms": bsms, "height": spec.real.max(axis=1)})
    for i, rng in enumerate(ppm_ranges):
        area = nt.integrate(ppm, spec, rng, part="complex")
        table[f"area_re_{i}"] = area.real
        table[f"area_im_{i}"] = area.imag
        table[f"center_{i}"] = nt.mean_ppm(ppm, spec, rng)
        table[f"width_{i}"] = nt.rms_width(ppm, spec, rng)
    return ppm, spec, table


def field_axes(table, bsms_zero_field=None):
    """
    Reference bsms values and ΔB columns, as in the original script.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`analyse_kpopt`.
    bsms_zero_field : float, optional
        bsms taken as ΔB = 0; the bsms of the tallest spectrum if ``None``.

    Returns
    -------
    table : pandas.DataFrame
        Copy of ``table`` with the original-script columns in front:
        ``integrals``, ``center``, ``width`` (first range) and the centred
        bsms / ΔB columns.
    summary : dict
        ``bsmsofmax``, ``bsmsofmincenter``, ``bsmsofmaxheight``,
        ``mT_per_bsms``.
    central_bsms : float
        bsms used as ΔB = 0.
    """
    last = len([c for c in table if c.startswith("center_")]) - 1
    summary = {
        "bsmsofmax": table.bsms[table.area_re_0.idxmax()],
        # The original script used the *last* ppm range here; kept as is.
        "bsmsofmincenter": table.bsms[table[f"center_{last}"].idxmin()],
        "bsmsofmaxheight": table.bsms[table.height.idxmax()],
        "mT_per_bsms": MT_PER_BSMS,
    }
    central = summary["bsmsofmaxheight"] if bsms_zero_field is None else bsms_zero_field
    legacy = pd.DataFrame({
        "integrals": table.area_re_0,
        "bsms": table.bsms,
        "center": table.center_0,
        "width": table.width_0,
        "bsms_centred_area": table.bsms - summary["bsmsofmax"],
        "bsms_centred_area_mT": (table.bsms - summary["bsmsofmax"]) * MT_PER_BSMS,
        "bsms_centred": table.bsms - central,
        "bsms_centred_mT": (table.bsms - central) * MT_PER_BSMS,
    })
    return pd.concat([legacy, table.drop(columns="bsms")], axis=1), summary, central


#%% ========================== FIGURES =====================================
def plot_diagnostics(table, ppm_ranges, title):
    """
    Integral (real, imaginary, modulus), centre and width of every region vs bsms.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`field_axes`.
    ppm_ranges : list of tuple
        The regions, for the column titles.
    title : str
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    n = len(ppm_ranges)
    fig, axs = plt.subplots(3, n, figsize=(3 * n, 7), sharex=True, squeeze=False)
    for i, rng in enumerate(ppm_ranges):
        re, im = table[f"area_re_{i}"], table[f"area_im_{i}"]
        axs[0, i].plot(table.bsms, re, "bo-", ms=3, label="real")
        axs[0, i].plot(table.bsms, im, "ro-", ms=3, label="imag")
        axs[0, i].plot(table.bsms, np.hypot(re, im), "ko-", ms=3, label="abs")
        axs[1, i].plot(table.bsms, table[f"center_{i}"], "o-", ms=3)
        axs[2, i].plot(table.bsms, table[f"width_{i}"], "o-", ms=3)
        axs[0, i].set_title(f"{rng[0]} to {rng[1]} ppm")
        axs[2, i].set_xlabel("bsms")
    axs[0, 0].set_ylabel("integral")
    axs[1, 0].set_ylabel("centre (ppm)")
    axs[2, 0].set_ylabel("RMS width (ppm)")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_summary(ppm, spec, table, ppm_ranges, title):
    """
    Integrals vs ΔB next to the spectral map (ppm vs ΔB).

    Parameters
    ----------
    ppm : numpy.ndarray
        Calibrated ppm axis.
    spec : numpy.ndarray
        Field-corrected spectra, one row per bsms.
    table : pandas.DataFrame
        Output of :func:`field_axes`.
    ppm_ranges : list of tuple
        The regions, for the legend.
    title : str
        Map title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax_int, ax_map : matplotlib.axes.Axes
    """
    dB = table.bsms_centred_mT.to_numpy()
    order = np.argsort(dB)
    fig, (ax_int, ax_map) = plt.subplots(1, 2, figsize=(10, 5), sharey=True,
                                         gridspec_kw={"width_ratios": [1, 3]})
    for i, marker in SUMMARY_RANGES:
        lo, hi = ppm_ranges[i]
        ax_int.plot(table[f"area_re_{i}"], dB, f"{marker}-", color="k",
                    label=f"integral {lo} to {hi} ppm")
    ax_int.set(xlabel="Integral [a.u.]", ylabel=r"$\Delta B$ (mT)")
    ax_int.grid(alpha=.3)
    ax_int.legend(fontsize=8, loc="upper right")

    real = spec.real[order]
    pcm = ax_map.pcolormesh(ppm, dB[order], real, shading="auto", cmap="inferno",
                            vmax=real.max() * MAP_VMAX_FRACTION)
    ax_map.set(xlim=(max(PLOT_RANGE), min(PLOT_RANGE)), ylim=(dB.max(), dB.min()),
               xlabel=r"$\delta$ (ppm)", title=title)
    ax_map.set_ylabel(r"$\Delta B$ (mT)")
    ax_map.tick_params(labelbottom=True, labelleft=True)
    ax_map.axhline(0, color="white", alpha=0.75, ls="--", lw=0.5)
    fig.colorbar(pcm, ax=ax_map, label="NMR signal [a.u.]")
    fig.tight_layout()
    return fig, ax_int, ax_map


def plot_slices(ppm, spec, table, fields_mT, mark_axes=()):
    """
    Spectra at the rows closest to the requested ΔB values.

    Parameters
    ----------
    ppm : numpy.ndarray
        Calibrated ppm axis.
    spec : numpy.ndarray
        Field-corrected spectra.
    table : pandas.DataFrame
        Output of :func:`field_axes`.
    fields_mT : list of float
        Requested ΔB values in mT.
    mark_axes : sequence of matplotlib.axes.Axes, optional
        Axes (e.g. of :func:`plot_summary`) where a dashed line of the same
        colour marks each plotted field.

    Returns
    -------
    matplotlib.figure.Figure
    """
    dB = table.bsms_centred_mT.to_numpy()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(6, 3))
    for idx, target in enumerate(fields_mT):
        k = int(np.argmin(np.abs(dB - target)))
        color = colors[idx % len(colors)]
        ax.plot(ppm, spec[k].real, lw=2, color=color, label=f"$\\Delta B$ = {dB[k]:.3g} mT")
        for other in mark_axes:
            other.axhline(dB[k], color=color, ls="--", alpha=.8)
    ax.set_xlim(max(PLOT_RANGE), min(PLOT_RANGE))
    ax.set(xlabel=r"$\delta$ (ppm)", ylabel="NMR signal [a.u.]")
    ax.grid(alpha=.3)
    for x in SLICE_GUIDES_PPM:
        ax.axvline(x, ls="--", lw=1, color="grey")
    fig.tight_layout()
    return fig


#%% ========================== OUTPUT ======================================
def save_outputs(savepath, ppm, spec, table, summary):
    """
    Write the CSV tables and one ``.dat`` file per spectrum.

    Parameters
    ----------
    savepath : str
        Output folder (created if needed); spectra go to ``spectra/``.
    ppm, spec : numpy.ndarray
        Calibrated axis and field-corrected spectra.
    table : pandas.DataFrame
        Output of :func:`field_axes`.
    summary : dict
        Output of :func:`field_axes`.
    """
    os.makedirs(os.path.join(savepath, "spectra"), exist_ok=True)
    table.to_csv(os.path.join(savepath, "integrals.csv"), index=False)
    pd.DataFrame([summary]).to_csv(os.path.join(savepath, "integrals_results.csv"), index=False)
    for k in range(spec.shape[0]):
        header = f"# BSMS: {table.bsms_centred[k]}, DeltaB: {table.bsms_centred_mT[k]}"
        np.savetxt(os.path.join(savepath, "spectra", f"spec{k:03d}.dat"),
                   np.column_stack([ppm, spec[k].real, spec[k].imag]), header=header)


#%% ========================== RUN =========================================
#: Results per experiment label: ppm, spec, table, summary, central_bsms.
results = {}

for exp in EXPERIMENTS:
    expdir = os.path.join(PATH, str(exp["expno"]))
    ds = nt.load(expdir)
    bsms = read_bsms_field(expdir)

    ppm, spec, table = analyse_kpopt(ds, bsms, PPM_RANGES, exp["ppm_offset"])
    table, summary, central = field_axes(table, exp["bsms_zero_field"])
    results[exp["label"]] = dict(ppm=ppm, spec=spec, table=table,
                                 summary=summary, central_bsms=central)

    title = f"cell: {exp['label']}; expn {exp['expno']}"
    plot_diagnostics(table, PPM_RANGES, title)
    fig_map, ax_int, ax_map = plot_summary(ppm, spec, table, PPM_RANGES, title)

    savepath = os.path.join(SAVEPATH_BASE, exp["label"])
    if SAVE:
        save_outputs(savepath, ppm, spec, table, summary)
        fig_map.savefig(os.path.join(savepath, "2dMap_figure.png"))

    ax_map.set_ylim(*ZOOM_MT)
    ax_int.set_ylim(*ZOOM_MT)
    fig_slices = plot_slices(ppm, spec, table, FIELDS_TO_PLOT_MT, (ax_int, ax_map))
    if SAVE:
        fig_map.savefig(os.path.join(savepath, "2dMap_figure_zoom.png"))
        fig_slices.savefig(os.path.join(savepath, "slices.png"))

plt.show()
