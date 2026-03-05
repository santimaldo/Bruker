# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate
import re
min_contact_time = 0 # discard p15<min_contact_time 

# # directorio de datos
# expns = np.concatenate([np.arange(60, 65), np.arange(66, 71)])
# path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-02_LiTFSI1M-aq_CA-cycles/"
# muestra = ""
# save = False
# plotRange = [100, -100]
# # rango de integracion
# ppmRanges = [[15, -10]
#             #[300, 150],
#             #[-0.5, -9]            
#             ]

# directorio de datos
expns = np.arange(20, 36)
path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2026-01-21_3.2mm_LiH/"
min_contact_time = 0
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\Debashis\Analysis\2026-01_LiH\CPspec/"
muestra = ""
save = False
plotRange = [100, -100]
# rango de integracion
ppmRanges = [[18, -10]
            #[300, 150],
            #[-0.5, -9]            
            ]

#=====================================================================
# 2D experiments
#=====================================================================
contact_times = np.zeros(len(expns))
Signals = np.zeros(len(expns))
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)

for jj, expn in enumerate(expns):
    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados(f'{path}/{expn}/',
                              read_pp = False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    contact_time = datos.acqus.dic["P"][15]
    if contact_time<min_contact_time:
        continue
    contact_times[jj] = contact_time
    ax_spec.plot(ppmAxis, spec)

    ###### start integrating the spectra
    colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
    ii = -1
    for ppmRange in ppmRanges:
        ii += 1
        color = colors[ii]
        ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
        r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
        ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
        ax_spec.axhline(0, color='k')

        signal = datos.Integrar(ppmRange=ppmRange)
        #tau_fit, signal_fit, residuals = datos.T1fit()
        Signals[jj] = signal

fig_popt, ax_popt = plt.subplots(num=382910)
ax_popt.plot(contact_times, Signals, 'o')#, color=color, label="Rising Edge")
#%%   
# guardo data:
if save:
    CPdata = np.array([contact_times, Signals]).T
    CPdata = CPdata[CPdata[:,1]!=0]
    CPdata = CPdata[CPdata[:,0].argsort()]

    header = "contact_time [s]\t"
    for ppmRange in ppmRanges:
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{muestra}_CP-curve.dat", CPdata, header=header)

# %%
