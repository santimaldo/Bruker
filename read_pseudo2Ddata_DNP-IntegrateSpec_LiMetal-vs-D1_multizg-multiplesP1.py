# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate
import matplotlib.cm as cm
import matplotlib.colors as mcolors


#=====================================================================
# Experimentos
#=====================================================================

expns, sample = np.arange(40, 63), 'Li dendrites + KBr'
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phases0_vsP1/"
Ns = np.ones_like(expns) * 2

# numero de puntos (RDs)
npts = 36

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-12-30_3.2mm_LiMetal-NO-DNP/"

# parametros de integracion espectral
plotRange = [350, 320]


#=====================================================================
# Loop principal
#=====================================================================

AMPS = []
p1s = []

for jj, expn in enumerate(expns):

    datos = DatosProcesados2D(f'{path}/{expn}/')

    # parametros del experimento
    p1s.append(datos.acqus.P1)

    datos.vdlist_path = f"lists/vd/sm2974.D1"
    recycle_delay = datos.get_vdlist(units='s')
    recycle_delay = recycle_delay[:npts]

    # espectro 1D
    datos.espectro.ppmSelect(plotRange)
    ppm = datos.espectro.ppmAxis
    re  = datos.espectro.real
    im  = datos.espectro.imag

    # integral compleja del espectro
    spec_cplx = re + 1j * im
    signal = integrate.simpson(spec_cplx, x=ppm, axis=1)
    signal = signal[:npts]

    AMPS.append(signal)

    # diagnostico
    if jj == 1:
        fig, ax = plt.subplots()
        ax.plot(ppm, re[jj] / np.max(re[jj]), 'k', lw=2)
        ax.set_xlim(np.max(ppm), np.min(ppm))
        ax.set_xlabel("NMR shift [ppm]")
        ax.set_ylabel("Norm. intensity")
        ax.set_title(f"expn {expn}")


#=====================================================================
# Ordenar por P1 creciente
#=====================================================================

p1s = np.array(p1s)
AMPS = np.array(AMPS)

idx = np.argsort(p1s)
p1s_sorted = p1s[idx]
AMPS_sorted = AMPS[idx]

ncurves = len(p1s_sorted)


#=====================================================================
# Colormap DISCRETO y UNIFORME
#=====================================================================

# colores igualmente espaciados en viridis
colors = cm.viridis(np.linspace(0, 1, ncurves))
colors = cm.turbo(np.linspace(0, 1, ncurves))

cmap = mcolors.ListedColormap(colors)

# indices enteros -> una celda por curva
bounds = np.arange(ncurves + 1) - 0.5
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# posiciones de los ticks (centros)
tick_locs = np.arange(ncurves)


#=====================================================================
# Plot: amplitud
#=====================================================================

fig_amps, ax_amps = plt.subplots()

norm_amp = np.abs(AMPS_sorted[0][-1])

for ii, signal in enumerate(AMPS_sorted):
    ax_amps.plot(recycle_delay,
                 np.abs(signal) / norm_amp,
                 'o-',
                 color=colors[ii])

ax_amps.set_xlabel("Recycle Delay [s]")
ax_amps.set_ylabel("Signal amplitude [a.u.]")
ax_amps.set_xscale('log')
ax_amps.set_yscale('log')
ax_amps.grid(True)

# colorbar uniforme
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig_amps.colorbar(sm,
                          ax=ax_amps,
                          ticks=tick_locs,
                          pad=0.02)

cbar.set_label(r"$P_1$ [$\mu$s]")
cbar.ax.set_yticklabels([f"{p1:.3f}" for p1 in p1s_sorted])

fig_amps.tight_layout()


#%%=====================================================================
# Plot: fase
#=====================================================================

fig_phase, ax_phase = plt.subplots()
zero_phase = np.mean( (np.angle(AMPS_sorted[:,-1])* 180 / np.pi) % 360 )
for ii, signal in enumerate(AMPS_sorted):
    phase = np.angle(signal) * 180 / np.pi
    phase = phase%360 - zero_phase
    phase = (phase + 1)%360 - 1
    print("WARNING! HARDCODED PHASE CORRECTION")
    ax_phase.plot(recycle_delay,
                  phase,
                  'o-',
                  color=colors[ii])

ax_phase.set_xlabel("Recycle Delay [s]")
ax_phase.set_ylabel(r"Phase [$^\circ$]")
ax_phase.set_xscale('log')
ax_phase.grid(True)

# colorbar uniforme
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig_phase.colorbar(sm,
                          ax=ax_phase,
                          ticks=tick_locs,
                          pad=0.02)

cbar.set_label(r"$P_1$ [$\mu$s]")
cbar.ax.set_yticklabels([f"{p1:.3f}" for p1 in p1s_sorted])

fig_phase.tight_layout()


# %%
