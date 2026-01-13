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
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import VoigtFit as vf




# metal
# expns = np.concatenate([np.arange(1, 60), np.arange(61, 100)])  # directorio de datos
expns, sample = np.arange(20,29), 'Li dendrites + KBr'
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phase90then0/"
Ns = np.arange(2,11)
labels = [f'N = {N}' for N in Ns]

expns, sample = np.arange(30,39), 'Li dendrites + KBr'
expns[0] = 20
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phase90_then-loop-of-0-180/"
Ns = np.arange(2,11)
labels = [f'N = {N}' for N in Ns]

expns, sample = np.arange(40,63), 'Li dendrites + KBr'
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phases0_vsP1/"
Ns = np.ones_like(expns) * 2

npts = 38 # puntos a promediar
nmeans = 20

absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-12-30_3.2mm_LiMetal-NO-DNP/"
# directorio de guradado
muestra = ""


#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
AMPS = []
p1s = []
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.set_fid()
    p1s.append(datos.acqus.P1)
    # extraigo los primeros npts de la fid
    fid0 = np.mean(datos.fid.real[:,:nmeans], axis=1) + 1j*np.mean(datos.fid.imag[:, :nmeans], axis=1)
    fid0 = fid0[:npts]

    if jj==1:
        slice = 30
        fig0, ax0 = plt.subplots()
        t = datos.acqus.DW * np.arange(datos.fid.real[0,:].size)
        ax0.plot(t*1000, datos.fid.real[slice,:], 'o-', label='Real')
        ax0.plot(t*1000, datos.fid.imag[slice,:], 'o-', label='Imaginary')
        ax0.axvspan(t[0]*1000, t[nmeans]*1000, color='gray', alpha=0.5, label="averaged")
        # ax0.set_xlim(-t[1]*1000, t[4*nmeans]*1000)
        ax0.set_xlabel("Acquisition time [ms]")
        ax0.legend()

    datos.vdlist_path = f"lists/vd/sm2974.D1"
    recycle_delay = datos.get_vdlist(units='s')
    recycle_delay = recycle_delay[:npts]
    AMPS.append(fid0)
    

    title = f"expn: {expn}\n"
    #ax_amps.plot(recycle_delay, np.abs(fid0), 'o-', label=labels[jj])
    

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
zero_phase = np.mean( (np.angle(AMPS_sorted[:,-1])* 180 / np.pi)  )
for ii, signal in enumerate(AMPS_sorted):
    phase = np.angle(signal) * 180 / np.pi
    phase = (phase - zero_phase+50)%360-50

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
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig_phase.colorbar(sm,
                          ax=ax_phase,
                          ticks=tick_locs,
                          pad=0.02)

cbar.set_label(r"$P_1$ [$\mu$s]")
cbar.ax.set_yticklabels([f"{p1:.3f}" for p1 in p1s_sorted])

fig_phase.tight_layout()


# %%
