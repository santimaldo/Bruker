# -*- coding: utf-8 -*-
"""
Created on Wed 12/3/2025

@author: Santi
"""
# raise Warning("NO FINCIONa!!!!!!!!!! CORREGIR")
from paths import NMRDATA_LOCAL, PROJECTS
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from Datos import *

# Data paths
expn = 19888
path = rf"{NMRDATA_LOCAL}\400dnp\2026-08-26_1.3mm_Gabriel_gf410_46\{expn}/"
savepath = rf"{PROJECTS}\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir"
muestra = "7Li_supercap_LiTFSI1M-aq_-1V"
save = False
plotRange = [280, -30]
Nsegments = 100


datos = DatosProcesadosT1(path, vdlist_path="lists/vd/sm2974")
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

# obtengo el espectro para el ultimo variable delay:
re = datos.espectro.real[-1]
im = datos.espectro.imag[-1]

fig_spec, ax_spec = plt.subplots(num=17856)
ax_spec.plot(ppmAxis, re)
ax_spec.set_xlim(np.max(plotRange), np.min(plotRange))


tau, signal = datos.get_T1data(plotRange)
tau = tau/1000


#%%============================================================================
# Lets try the pseudo 2d exponential fitting
#=============================================================================

Nind, Ndir = spec.shape
ydata2d = spec

# Compute the number of columns per segment
columns_per_segment = Ndir // Nsegments  # Integer division
# Trim extra columns if needed
trimmed_columns = columns_per_segment * Nsegments
ydata2d_trimmed = ydata2d[:, :trimmed_columns]  # Remove excess columns
# Reshape and compute the mean over the last dimension
ydata2d_reduced = ydata2d_trimmed.reshape(Nind, Nsegments, columns_per_segment).mean(axis=2)
ppmAxis_reduced = np.linspace(plotRange[0], plotRange[1], Nsegments)


# funcion de decaimiento
def T1fit(t, y0, A, T1):
    return y0 + A * np.exp(-t/T1)


T1_values = np.zeros(Nsegments)
T1_errors = np.zeros(Nsegments)
R2_values = np.zeros(Nsegments)


for ii in range(Nsegments):
    print(f"Segment {ii+1} of {Nsegments}")

    ydata = ydata2d_reduced[:, ii]

    # Initial guess
    y0_guess = ydata[-1]
    A_guess = ydata[0] - y0_guess
    T1_guess = (tau[-1] - tau[0]) / 2

    try:
        popt, pcov = curve_fit(
            T1fit, tau, ydata,
            p0=[y0_guess, A_guess, T1_guess],
            bounds=([-np.inf, -np.inf, 0], [np.inf, np.inf, np.inf])
        )

        y0, A, T1 = popt

        T1_values[ii] = T1
        T1_errors[ii] = np.sqrt(np.diag(pcov))[2]

        # R2 del ajuste
        residuals = ydata - T1fit(tau, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata - np.mean(ydata))**2)

        if ss_tot > 0:
            R2_values[ii] = 1 - ss_res/ss_tot
        else:
            R2_values[ii] = np.nan

    except (RuntimeError, ValueError):
        T1_values[ii] = np.nan
        T1_errors[ii] = np.nan
        R2_values[ii] = np.nan


#%%

fig, ax = plt.subplots()

ax.pcolormesh(
    ppmAxis_reduced,
    tau,
    ydata2d_reduced
)

ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
ax.set_xlabel("ppm")
ax.set_ylabel("tau [s]")


#%%

fig, axs = plt.subplots(ncols=2, figsize=(10, 5))

ax = axs[0]
ax.plot(ydata2d_reduced.T, 'o-')
ax.set_xlabel("ppm axis")
ax.set_ylabel("Signal")

ax = axs[1]
ax.plot(tau, ydata2d_reduced, 'o-')
ax.set_xlabel("tau axis")
ax.set_ylabel("Signal")


#%% PPM-T1 map with projections
#########################################################
#########################################################
#########################################################
#########################################################

from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(10, 7))

gs = GridSpec(
    2, 2,
    width_ratios=[4, 1],
    height_ratios=[1, 4],
    hspace=0.05,
    wspace=0.05
)

# -------------------------------------------------------------------------
# Axes
# -------------------------------------------------------------------------
ax_top = fig.add_subplot(gs[0, 0])
ax_map = fig.add_subplot(gs[1, 0], sharex=ax_top)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax_map)

# -------------------------------------------------------------------------
# Main signal map
# -------------------------------------------------------------------------
pcm = ax_map.pcolormesh(
    ppmAxis_reduced,
    tau,
    ydata2d_reduced
)

ax_map.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))

ax_map.set_xlabel("ppm")
ax_map.set_ylabel("tau [s]")

# -------------------------------------------------------------------------
# Vertical projection -> TOP
# Integrate/sum over tau
# -------------------------------------------------------------------------
projection_ppm = np.sum(ydata2d_reduced, axis=0)

ax_top.plot(
    ppmAxis_reduced,
    projection_ppm
)

ax_top.set_xlim(
    np.max(ppmAxis_reduced),
    np.min(ppmAxis_reduced)
)

ax_top.set_ylabel("Projection")

# Hide x tick labels because they are shown on the main map
ax_top.tick_params(axis="x", labelbottom=False)

# -------------------------------------------------------------------------
# Horizontal projection -> RIGHT
# Integrate/sum over ppm
# -------------------------------------------------------------------------
projection_tau = np.sum(ydata2d_reduced, axis=1)

ax_right.plot(
    projection_tau,
    tau
)

# Hide y tick labels because they are shown on the main map
ax_right.tick_params(axis="y", labelleft=False)

ax_right.set_xlabel("Projection")

# -------------------------------------------------------------------------
# Optional: make projection axes start at zero
# -------------------------------------------------------------------------
ax_top.set_ylim(bottom=0)
ax_right.set_xlim(left=0)


#%% T1 vs ppm
#########################################################
#########################################################
#########################################################

fig, ax = plt.subplots(figsize=(10, 5))

ax.errorbar(
    ppmAxis_reduced,
    T1_values,
    yerr=T1_errors,
    fmt='o-'
)

ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
ax.set_xlabel("ppm")
ax.set_ylabel("T1 [s]")

plt.show()


#%% T1 vs ppm without error bars
#########################################################
#########################################################

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(
    ppmAxis_reduced,
    T1_values,
    'o-'
)

ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
ax.set_xlabel("ppm")
ax.set_ylabel("T1 [s]")
ax.set_ylim(0, 1)

plt.show()


#%% R2 vs ppm
#########################################################
#########################################################

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(
    ppmAxis_reduced,
    R2_values,
    'o-'
)

ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
ax.set_xlabel("ppm")
ax.set_ylabel("R²")
ax.set_ylim(0, 1)

plt.show()


#%%

# # guardo data:
# if save:
#     for ii, fig in enumerate(figures):
#         filename = f'{savepath}/{muestra}_T1_range{ii}.png'
#         fig.savefig(filename)   # save the figure to file

#     Signals = np.array(Signals).T
#     tau = tau.reshape(tau.size, 1)
#     T1data = np.hstack((tau, Signals))
#     header = "tau [s]\t"
#     for ii, ppmRange in enumerate(ppmRanges):
#         header += f"{ppmRange} ppm\t"
#     np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

#     data = np.array([ppmAxis, re, im]).T
#     np.savetxt(f"{savepath}/{muestra}_lastSpectrum.dat", data)
#     fig_spec.savefig(f"{savepath}/{muestra}_lastSpectrum_withRanges.png")


#%% T1 vs ppm without error bars
#########################################################
#########################################################

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(
    ppmAxis_reduced,
    T1_values*(ppmAxis_reduced)**2,
    'o-'
)

ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
ax.set_xlabel("ppm")
ax.set_ylabel(r"$T_1 K^2$ [a.u]")
# ax.set_ylim(0, 1)

plt.show()
# %%
