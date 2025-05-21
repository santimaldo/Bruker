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

# directorio de datos
expns = [4,5]#,6]
initial_times_h = [0, 29.75]#, 42.28]
expn_durations_h = [29.47, 12.43]#, 0.56667]

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-04-29_IMEC_eLi-sym_LP40-OCV/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_in-situ_LP40\00_OCV/"
muestra = "7Li_eLi-symm-cell_LP40"

save = False
plotRange = [300, 200]
# rango de integracion
ppmRange = [300,200]          
            

#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
# datos = DatosProcesados(f'{path}/{expn_before}/')
# datos.espectro.ppmSelect(plotRange)
# ppmAxis = datos.espectro.ppmAxis
# re = datos.espectro.real
# ppmAxis = datos.espectro.ppmAxis
# # find the ppm of the maximum in the range < ppm_trershold ppm
# re_in_ROI = re[ppmAxis < ppm_treshold] # re in Region Of Interest
# ppmAxis_in_ROI = ppmAxis[ppmAxis < ppm_treshold]
# max_index = np.argmax(re_in_ROI)
# ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]

colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)
fig,ax = plt.subplots(num=1785731, figsize=(6, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(6, 3))
fig_int.suptitle(muestra)
for jj, expn in enumerate(expns):
    #=====================================================================
    # Ajuste de espectro 2D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    # Remove rows in the spectrum where there are only zeros
    spec = spec[np.sum(spec, axis=1)!=0, :]

    tau = np.linspace(0, expn_durations_h[jj], spec.shape[0]) + initial_times_h[jj]  # horas

    # grafico todos los espectros juntos
    ppm_of_max = []
    for ii in range(datos.espectro.size[0]):
        re = datos.espectro.real[ii]
        im = datos.espectro.imag[ii]
        if expn==4 or expn==6:
            ax_spec.plot(ppmAxis, re,
                        color=colors[jj], alpha=0.2)# 0.2+(ii/datos.espectro.size[0])*0.7)
        # find the ppm of the maximum in the range < ppm_treshold ppm
        re_in_ROI = re # re in Region Of Interest
        ppmAxis_in_ROI = ppmAxis
        max_index = np.argmax(re_in_ROI)
        ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]
        ppm_of_max.append(ppm_of_max_in_equilibrium)
    ppm_of_max = np.array(ppm_of_max)
    ppm_of_max = ppm_of_max[:tau.size]

    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvline(r1, color='k', linestyle='--')
    ax_spec.axvline(r2, color='k', linestyle='--')
    ax_spec.set_xlabel('chemical shift [ppm]')
    ax_spec.set_ylabel('Intensity [a.u.]')
    ax_spec.axhline(0, color='k')


    ax.plot(tau, ppm_of_max, 'o-', color=colors[0])


    signal = datos.Integrar(ppmRange=ppmRange)
    signal = signal[:tau.size]
    if jj == 0:
        initial_signals = signal[0]
    signal = signal/initial_signals  # normalizo la señal al primer espectro

    # -------------
    ax_int.plot(tau, signal[:tau.size], 'o', color=colors[0])
    if expn==5:
        # firts plating
        fpd = 0.8 # first plating duration in hours
        condition = tau<tau[0]+fpd
        tau_plating = tau[condition]

        signal_plating = signal[condition]
        ax_int.plot(tau_plating, signal_plating, 'o', color=colors[jj])
        ppm_of_max_plating = ppm_of_max[condition]
        ax.plot(tau_plating, ppm_of_max_plating, 'o-', color=colors[jj])
        for ii in range(tau.size):
            if tau[ii] < tau[0]+fpd:
                color = colors[jj] 
            else:
                color = colors[0]   
            re = datos.espectro.real[ii]
            ax_spec.plot(ppmAxis, re,   
                        color=color, alpha=0.2)
        # second plating
        ocptbp = 8 # ocp time between the two plating
        spd = 0.4 # second plating duration in hours
        condition = (tau>tau[0]+fpd+ocptbp) &(tau<tau[0]+fpd+ocptbp+spd)
        tau_plating = tau[condition]
        signal_plating = signal[condition]
        ax_int.plot(tau_plating, signal_plating, 'o', color=colors[jj])
        ppm_of_max_plating = ppm_of_max[condition]
        ax.plot(tau_plating, ppm_of_max_plating, 'o-', color=colors[jj])
        for ii in range(tau.size):
            if (tau[ii]>tau[0]+fpd+ocptbp) & (tau[ii]<tau[0]+fpd+ocptbp+spd):    
                re = datos.espectro.real[ii]
                ax_spec.plot(ppmAxis, re,   
                            color=colors[jj], alpha=0.15)



ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')
ax_spec.set_xlim(280,230)

ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta_{max}$ [ppm]')
ax.set_ylim(254.5, 256.5)




# guardo data:
if save:
    filename = f'{savepath}/{muestra}_T1.png'
    fig.savefig(filename)   # save the figure to file

    Signals = np.array(Signals).T
    tau = tau.reshape(tau.size, 1)
    T1data = np.hstack((tau, Signals))
    header = "tau [s]\t"
    for ppmRange in ppmRanges:
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)
