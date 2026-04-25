# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi



CORREGIR: TIEMPOS Y FASES
"""

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate

########### R3 - PC
# directorio de datos
expns = np.arange(10, 300)
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2026-04-22_safebatt_Gr-NMC_cellG3/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\SafeBatt\GraphiteNMC\Analysis\2026-04_cellG3_NMR/"
muestra = "Cell_G3"
save = True
# rango de guardado
ppmRange = [600, -200]


################################################
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))
tau = np.zeros(expns.size)

for jj, expn in enumerate(expns):
    expn = int(expn)    
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    ppmAxis = datos.espectro.ppmAxis
    re = datos.espectro.real
    im = datos.espectro.imag

    spec_time = datos.acqus.dic['DATE_START']
    # spec_time = datos.acqus.dic['DATE_START'] - Nominal_duration_of_experiment # to take into account the time of autotune
    if jj==0:
        t_0 = spec_time# s
        spec_vs_t = np.zeros([expns.size, re.size])
    tau[jj] = (spec_time - t_0)/3600  # convert to hours
    

    spec_vs_t[jj, :] = re

    ax_spec.plot(ppmAxis, re)# 0.2+(ii/datos.espectro.size[0])*0.7)
    
    np.savetxt(f'{savepath}/1dspectra/spec_{jj:03d}.dat',
               np.array([ppmAxis, re]).T,
               header='ppm, Intensity [a.u.]')

#%% ax_spec.

np.savetxt(f'{savepath}/1dspectra/expn_list.dat',
            expns,
            header='list of experiment numbers')
np.savetxt(f'{savepath}/1dspectra/time_list.dat',
            tau,
            header='time [h]')

#%%
fig, ax = plt.subplots()
zlim = [0, 2e6]
ax.pcolormesh(ppmAxis,  tau, spec_vs_t, vmin=min(zlim), vmax=max(zlim))
ax.set_xlabel("ppm")
ax.set_ylabel("Time [h]")
# ax.set_xlim(max(ppmRange), min(ppmRange))
ax.set_xlim([100, -100])
# %%
