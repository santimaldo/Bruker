# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi



CORREGIR: TIEMPOS Y FASES
"""

import os

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate

# ########### G3 - NMC-Graphite
# # directorio de datos
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-04-22_safebatt_Gr-NMC_cellG3/"
# # zgs
# files_expns = [[10,335,1],
#                [337,565,2]]  
# expns = np.array([])
# nominal_durations = np.array([])
# nominal_duration = [13*60+43, 6*60+53]  # seconds
# for j, [start, stop, step] in enumerate(files_expns):
#     expns_j = np.arange(start, stop+1, step)
#     expns = np.append(expns, expns_j)
#     nominal_durations = np.append(nominal_durations, np.ones_like(expns_j)*nominal_duration[j]) # s, to take into account the time of autotune
# # # # hahnechos
# # # files_expns = [336,564,2] 
# # # expns = np.arange(start, stop+1, step)
# # # expns = expns[expns != 444] # remove expn 444 which is an outlier
# # # nominal_durations = np.ones_like(expns)*(13*60+43) # s, to take into account the time of autotune
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\SafeBatt\Analysis\2026-04_cellG3_NMR/"
# muestra = "Cell_G3"#_hahnecho"
# save = True
# # rango de guardado
# ppmRange = [600, -200]



# ########### R2 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(10, 200)
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2025-08-10_ccATMC_Rui-R1_NMC-Cu_PC"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-08_R2/"
# muestra = "Cell_R2"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 14*60+18 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune



# # ########### R3 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(30, 330)
# expns = expns[expns%10 != 0] 
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2025-10-01_ccATMC_Rui-R3_NMC-Cu_PC"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-10_R3/"
# muestra = "Cell_R3"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 14*60+18 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune


# # ########### R5 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(10, 110)
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-04-05_ccATMC_Rui-R5_NMC-Cu_CC"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-04_in-situ_repeat__R5-R6\R5/"
# muestra = "Cell_R5"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 14*60+18 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune



# # ########### R6 - NMC-Cu CC protocol
# # directorio de datos
# expns = np.arange(10,71)
# path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-04-06_ccATMC_Rui-R6_NMC-Cu_CC"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-04_in-situ_repeat__R5-R6\R6/"
# muestra = "Cell_R6"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 14*60+18 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune



# # ########### R9 - NMC-Cu CC protocol
# # directorio de datos
# expnlist = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-05-07_ATMC1_Rui-R9_NMC-Cu_CC\0_datalists\expnlist.txt"
# expns = np.loadtxt(expnlist, dtype=int)
# path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-05-07_ATMC1_Rui-R9_NMC-Cu_CC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-05_R9_CC/"
# muestra = "Cell_R9"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 13*60+43 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune

# ########### R10 - NMC-Cu CC protocol
# # directorio de datos
# expnlist = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-05-13_ATMC1_Rui-R10_NMC-Cu_PC\0_datalists\expnlist.txt"
# expns = np.loadtxt(expnlist, dtype=int)
# path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-05-13_ATMC1_Rui-R10_NMC-Cu_PC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-05_R10_PC/"
# muestra = "Cell_R10"
# save = True
# # rango de guardado
# ppmRange = [149, -149]
# Nominal_duration = 13*60+43 # seconds
# nominal_durations = np.ones_like(expns)*Nominal_duration # s, to take into account the time of autotune


########### G6 - NMC-Gr Safebatt
# directorio de datos
litsts_path = r"C:\Users\Santi\Documents\NMRdata\300neo\2026-08-09_safebatt_Gr-NMC_cellG6\0_datalists/"
expns = np.loadtxt(litsts_path+"expnlist.txt", dtype=int)
path = rf"C:\Users\Santi\Documents\NMRdata\300neo\2026-08-09_safebatt_Gr-NMC_cellG6/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\SafeBatt\Analysis\2026-08_cellG6_NMR/"
folder_suffix = "_LiMetal_zg"
muestra = "Cell_G6"
save = True
# rango de guardado
ppmRange = [399, 201]
nominal_durations = np.loadtxt(litsts_path+"durations.txt")


################################################
# Crea la carpeta si no existe
os.makedirs(f"{savepath}/1dspectra{folder_suffix}", exist_ok=True)
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))
tau = np.zeros(expns.size)
time = np.zeros_like(tau)
for jj, expn in enumerate(expns):
    expn = int(expn)    
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    ppmAxis = datos.espectro.ppmAxis
    re = datos.espectro.real
    im = datos.espectro.imag

    
    # spec_time = datos.acqus.dic['DATE'] - Nominal_duration_of_experiment # to take into account the time of autotune
    Nominal_duration = nominal_durations[jj]
    spec_time       = datos.acqus.dic['DATE'] - Nominal_duration
    if jj==0:
        t_0 = spec_time # s
        spec_vs_t = np.zeros([expns.size, re.size])
    
    tau[jj] = (spec_time - t_0)
    
    spec_vs_t[jj, :] = re
    # speci_vs_t[jj, :] = im

    ax_spec.plot(ppmAxis, re)# 0.2+(ii/datos.espectro.size[0])*0.7)
    
    np.savetxt(f'{savepath}/1dspectra{folder_suffix}/spec_{jj:04d}.dat',
               np.array([ppmAxis, re, im]).T,
               header='ppm, Intensity [a.u.]')

#%% ax_spec.

np.savetxt(f'{savepath}/1dspectra{folder_suffix}/expn_list.dat',
            expns,
            header='list of experiment numbers')
np.savetxt(f'{savepath}/1dspectra{folder_suffix}/time_list.dat',
            tau,
            header='time [h]')

with open("info.txt", 'w') as f:    
    f.write(f"Nominal duration of experiments: {Nominal_duration} s = {Nominal_duration/60} min")

#%%
fig, ax = plt.subplots()
ax.pcolormesh(ppmAxis,  tau, spec_vs_t)#, vmin=min(zlim), vmax=max(zlim))
ax.set_xlabel("ppm")
ax.set_ylabel("Time [h]")
# ax.set_xlim(max(ppmRange), min(ppmRange))
# ax.set_xlim([100, -100])
# %%
