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



# # ########### R12 - NMC-Cu CC protocol
# # directorio de datos
# expnlist = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-28_ATMC1_Rui-R12_NMC-Cu_CC\0_datalists\expnlist.txt"
# expns = np.loadtxt(expnlist, dtype=int)
# path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-28_ATMC1_Rui-R12_NMC-Cu_CC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R12_CC/"
# muestra = "Cell_R12"
# folder_suffix = "_LiMetal_zg"
# save = True
# # rango de guardado
# ppmRange = [350, 200]


# # ########### R13 - NMC-Cu PC protocol
# # directorio de datos
# expnlist = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-31_ATMC1_Rui-R13_NMC-Cu_PC\0_datalists\expnlist.txt"
# expns = np.loadtxt(expnlist, dtype=int)
# path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-31_ATMC1_Rui-R13_NMC-Cu_PC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R13_PC/"
# muestra = "Cell_R13"
# folder_suffix = "_LiMetal_zg"
# save = True
# # rango de guardado
# ppmRange = [400, -200]


# ########### R14 - NMC-Cu CC protocol
# directorio de datos
expnlist = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-09-02_ATMC1_Rui-R14_NMC-Cu_CC\0_datalists\expnlist.txt"
expns = np.loadtxt(expnlist, dtype=int)
path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-09-02_ATMC1_Rui-R14_NMC-Cu_CC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-09_R14_CC/"
muestra = "Cell_R14"
folder_suffix = "_LiMetal_zg"
save = True
# rango de guardado
ppmRange = [350, 200]

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

    
    # Si asumo que todos los experimentos tienen la misma duracion, puedo usar directamente el tiempo en el que termino el experimento
    spec_time       = datos.acqus.dic['DATE']
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

np.savetxt(f'{savepath}/1dspectra{folder_suffix}/expn.list',
            expns,
            header='list of experiment numbers')
np.savetxt(f'{savepath}/1dspectra{folder_suffix}/time.list',
            tau,
            header='time [s]')
np.savetxt(f'{savepath}/1dspectra{folder_suffix}/initial_time.list',
            [t_0],
            header='time [s]')
#%%
fig, ax = plt.subplots()
ax.pcolormesh(ppmAxis,  tau, spec_vs_t)#, vmin=min(zlim), vmax=max(zlim))
ax.set_xlabel("ppm")
ax.set_ylabel("Time [h]")
# ax.set_xlim(max(ppmRange), min(ppmRange))
# ax.set_xlim([100, -100])
# %%
