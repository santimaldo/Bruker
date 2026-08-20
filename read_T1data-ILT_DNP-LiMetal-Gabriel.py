# -*- coding: utf-8 -*-
"""
Created on Wed 12/3/2025

@author: Santi
"""

import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Laplace import ILT

# Data paths
expn = 51
path = rf"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\2026-08-17_1.3mm_Gabriel_gf410_41\{expn}/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir"
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


# funcion de kernel
def T2(x, tt):
    return 1-np.exp(-x/tt)

Npts_L = 200
alpha = 1e-4

# Inicializo la clase ILT
# ilt = ILT(alpha=alpha, rango=(1e0, 1e4), kernel=T2, Nilt=Npts_L,
#             figure=1, savepath=savepath)

tau, signal = datos.get_T1data(plotRange)
tau = tau/1000

colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
Signals = []
ii = -1

figures = []
# for ppmRange in ppmRanges:
#     ii += 1
#     color = colors[ii]
#     ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
#     r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
#     ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
#     ax_spec.axhline(0, color='k')

#     tau, signal = datos.get_T1data(ppmRange)
#     signal_reordered = - (signal - signal[-1])
#     Signals.append(signal)

#     ilt.DoTheStuff(signal_reordered, tau, muestra=f"-1V_ppmRange_{ppmRange[0]:.2f}_to_{ppmRange[1]:.2f}")
#%%============================================================================
# Lets try the pseufo 2d ILTFT

# Inicializo la clase ILT
# ilt = ILT(alpha=alpha, rango=(0.05, 100), kernel=T2, Nilt=Npts_L,
            # figure=None, savepath=savepath)
# iltft, ydata2d_reduced = ilt.ILTFT(spec, tau, Nsegments=Nsegments)
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

Nilt = 100
iltft = np.zeros((Nilt, Nsegments))
for ii in range(Nsegments):
    print(f"Segment {ii+1} of {Nsegments}")
    ydata = ydata2d_reduced[:,ii]
    ilt = ILT(alpha=alpha, rango=(0.001, 1), kernel=T2, Nilt=Nilt, figure=ii, savepath=savepath)
    ilt.DoTheStuff(ydata, tau, muestra=f"{ppmAxis_reduced[ii]:.2f} ppm")
    iltft[:,ii] = ilt.yilt


#%%
fig, ax = plt.subplots()
ax.pcolormesh(ppmAxis_reduced, tau, ydata2d_reduced)
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


#%%
fig, ax = plt.subplots(figsize=(10, 5))
# ax.contour(ppmAxis_reduced, ilt.xilt, iltft, 5)#, vmax=0.05, vmin=0)
ax.pcolormesh(ppmAxis_reduced, ilt.xilt, iltft)#, vmax=0.1)
ax.set_yscale("log")
ax.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))
# ax.set_ylim(0.05, 1e2)
ax.set_xlabel("ppm")
ax.set_ylabel("T1 [s]")

#%%




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


#%% PPM-T2 map with projections
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
# Main ILTFT map
# -------------------------------------------------------------------------
pcm = ax_map.pcolormesh(
    ppmAxis_reduced,
    ilt.xilt,
    iltft,
    vmax=0.8 * iltft.max()
)

# ax_map.set_yscale("log")
ax_map.set_xlim(np.max(ppmAxis_reduced), np.min(ppmAxis_reduced))

ax_map.set_xlabel("ppm")
ax_map.set_ylabel("T1 [s]")
ax_map.set_ylim([0.05, 0.2])

# -------------------------------------------------------------------------
# Vertical projection -> TOP
# Integrate/sum over T2
# -------------------------------------------------------------------------
projection_ppm = np.sum(iltft, axis=0)

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
projection_T2 = np.sum(iltft, axis=1)

ax_right.plot(
    projection_T2,
    ilt.xilt
)

# ax_right.set_yscale("log")

# Hide y tick labels because they are shown on the main map
ax_right.tick_params(axis="y", labelleft=False)

ax_right.set_xlabel("Projection")

# -------------------------------------------------------------------------
# Optional: make projection axes start at zero
# -------------------------------------------------------------------------
ax_top.set_ylim(bottom=0)
ax_right.set_xlim(left=0)

plt.show()
# %%
