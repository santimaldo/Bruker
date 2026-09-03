# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi


"""

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate
import os



# # ############ R12 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(31, 1200, 10)
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-08-28_ATMC1_Rui-R12_NMC-Cu_CC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R12_CC/"
# muestra = "Cell_R12"
# save = True
# ppmRange = None #[2300, -1300]

# # ############ R13 - NMC-Cu PC protocol
# # directorio de datos
# expns = np.arange(31, 922, 10)
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-08-31_ATMC1_Rui-R13_NMC-Cu_PC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R13_PC/"
# muestra = "Cell_R13"
# save = True
# ppmRange = None #[2300, -1300]

# ############ R14 - NMC-Cu CC protocol
# directorio de datos
expns = np.arange(31, 351, 10)
path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-09-02_ATMC1_Rui-R14_NMC-Cu_CC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-09_R14_CC/"
muestra = "Cell_R14"
save = True
ppmRange = None #[2300, -1300]



################################################
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))

os.makedirs(f"{savepath}/VOCS", exist_ok=True)
tau = []
tau_zg = []
VOCS_vs_time = []
atmc_durations = []
for jj, expn in enumerate(expns):
    expn = int(expn)
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados2D(f'{path}/{expn}/', remove_empty = True)
    if ppmRange is not None:
        datos.espectro.ppmSelect(ppmRange)

    Nfreqs = datos.acqu2s.TD
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.spec

    #### Calculating time of the spectra:
    expn_time = datos.acqus.dic['DATE']
    
    # Loop over the spectra in the expn    
         
    VOCS = np.sum(np.abs(spec), axis=0) 
    
    tau.append([expn, expn_time])
    VOCS_vs_time.append(VOCS)

    np.savetxt(f'{savepath}/VOCS/VOCS_{expn:03d}.dat',
                np.array([ppmAxis, VOCS]).T,
                header='ppm\tIntensity[a.u.]')

# Saving time list
# Create output directories if they do not exist
np.savetxt(f'{savepath}/VOCS/time.list',
            np.array(tau),
            header='Experiment Number\tTime [s]')

# Convert lists to numpy arrays
VOCS_vs_time = np.array(VOCS_vs_time)

# Time axis [s]
times = np.array(tau)[:,1]

#%%
# =========================================================================
# Plots
# =========================================================================

# -------------------------------------------------------------------------
# Overlay of all VOCS spectra
# -------------------------------------------------------------------------
fig_spec, ax_spec = plt.subplots(num=2001, figsize=(6, 4))
maxx = 0.3e7
for spectrum in VOCS_vs_time:
    ax_spec.plot(ppmAxis, spectrum, lw=0.8)

ax_spec.set_xlabel(r"$^7$Li Chemical shift [ppm]")
ax_spec.set_ylabel("Intensity [a.u.]")
ax_spec.set_title("VOCS spectra")
ax_spec.set_ylim([0, maxx])
ax_spec.grid(alpha=0.3)

# -------------------------------------------------------------------------
# Heatmap: intensity vs ppm vs time
# -------------------------------------------------------------------------
ppm_mesh, time_mesh = np.meshgrid(ppmAxis, times)

fig_map, ax_map = plt.subplots(num=2002, figsize=(4, 6))

pcm = ax_map.pcolormesh(
    ppm_mesh,
    time_mesh,
    VOCS_vs_time,
    vmax=maxx,
    shading='nearest'
)

cbar = fig_map.colorbar(pcm, ax=ax_map)
cbar.set_label("Intensity [a.u.]")

ax_map.set_xlabel(r"$^7$Li Chemical shift [ppm]")
ax_map.set_ylabel("Time [s]")

ax_map.set_xlim([ppmAxis.max(), ppmAxis.min()])

plt.tight_layout()
# %%
