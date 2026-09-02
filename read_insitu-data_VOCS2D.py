# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi



ESTA HECHO PARA SOLO 1 ZG
"""

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate
import os



# ############ R9 - NMC-Cu CC protocol
# directorio de datos
expns = np.arange(1, 40)
path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-05-07_ATMC1_Rui-R9_NMC-Cu_CC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-05_R9_CC/"
muestra = "Cell_R9"
save = True
# rango de guardado
nominal_durations = {"VOCS": 166.07,  # seconds
                     "zg": 388.61}
spectra_per_expn = {"VOCS": 11,
                    "zg": 1,
                    }
ppmRange = None #[2300, -1300]
# ########### R10 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(1, 48)
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2026-05-13_ATMC1_Rui-R10_NMC-Cu_PC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-05_R10_PC/"
# muestra = "Cell_R10"
# save = True
# # rango de guardado
# nominal_durations = {"VOCS": 166.07,  # seconds
#                      "zg": 388.61}
# spectra_per_expn = {"VOCS": 11,
#                     "zg": 1,
#                     }
# ppmRange = None #[2300, -1300]



################################################
spectra_per_loop = spectra_per_expn["VOCS"] + spectra_per_expn["zg"]
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

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.spec

    #### Calculating time of the spectra:
    total_time = datos.acqus.dic['DATE'] - datos.acqus.dic['DATE_START']
    # Nspectra = datos.proc2s.SI
    # Nvocs = Nspectra // spectra_per_loop
    nominal_duration_total = nominal_durations["VOCS"] + nominal_durations["zg"]
    ATMC_time = total_time - nominal_duration_total # to take into account the time of autotune and the duration of the spectra
    initial_expn_time = datos.acqus.dic['DATE_START'] + ATMC_time
    atmc_durations.append(ATMC_time) 
    if jj==0:
        t_0 = initial_expn_time
    initial_expn_time = initial_expn_time - t_0 # time relative to the first spectrum
    
    # Loop over the spectra in the expn    
    for ii in range(datos.proc2s.SI//spectra_per_loop):        
        VOCS_ii = np.sum(np.abs(spec[ii*spectra_per_loop:ii*spectra_per_loop+spectra_per_expn['VOCS'], :]), axis=0) 
        time_VOCS_ii = initial_expn_time + ii*nominal_duration_total
        
        time_zg_ii = time_VOCS_ii + nominal_durations['VOCS']
        tau.append([expn, ii, time_VOCS_ii])
        tau_zg.append([int(f"{expn:02d}999{ii:03d}"),time_zg_ii])
        VOCS_vs_time.append(VOCS_ii)

        np.savetxt(f'{savepath}/VOCS/VOCS_{expn:02d}_{ii:02d}.dat',
                   np.array([ppmAxis, VOCS_ii]).T,
                   header='ppm\tIntensity[a.u.]')

# Saving time list
# Create output directories if they do not exist
np.savetxt(f'{savepath}/VOCS/time_list.dat',
            np.array(tau),
            header='Experiment Number\tSpectrum Number\tTime [h]')
np.savetxt(f'{savepath}/time_list_zg.dat',
            np.array(tau_zg),
            header='Experiment Number\tSpectrum Number\tTime [h]')


# Convert lists to numpy arrays
VOCS_vs_time = np.array(VOCS_vs_time)

# Time axis [s]
times = np.array([row[2] for row in tau])

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
    shading='auto'
)

cbar = fig_map.colorbar(pcm, ax=ax_map)
cbar.set_label("Intensity [a.u.]")

ax_map.set_xlabel(r"$^7$Li Chemical shift [ppm]")
ax_map.set_ylabel("Time [s]")

ax_map.set_xlim([ppmAxis.max(), ppmAxis.min()])

plt.tight_layout()
# %%
