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

# metal
expns = [1,2]# directorio de datos
plotRange = [500, 100]
ppmRange = [400,200]   # rango de integracion   

absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-06-17_3.2mm_IMECdendrites_KBr-T1-while-mw-ON/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\in-situ\2025-05_eLI-LFP_LP40/"
muestra = "eLi-KBr"


     
            

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


spec_vs_t_list = []
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
fig,ax = plt.subplots(num=1785731, figsize=(8, 3))
figph,axph = plt.subplots(num=178573111, figsize=(8, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(8, 3))
fig_int.suptitle(muestra)
for jj, expn in enumerate(expns):
    #=====================================================================
    # Ajuste de espectro 2D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)
    # True number of spectra:
    Nspec1d = datos.acqu2s.dic['TD']
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real[:Nspec1d, :]
    spec_vs_t = np.zeros((Nspec1d, spec.shape[1]))
    if jj==0:
        t_0 = datos.acqus.dic['DATE_START'] # s
        #### NOTE:
        # The time points should be shifted by DELTA, wich is the time between 1D slices, since it started with the waiting.
        # for simplicity, I will not do it here.
    t_ini = datos.acqus.dic['DATE_START']    
    t_end = datos.acqus.dic['DATE'] # s
    tau = np.linspace(t_ini, t_end, Nspec1d) - t_0 # horas
    tau = tau/3600  # convert to hours
    # grafico todos los espectros juntos
    ppm_of_max = []
    phases = []
    for ii in range(Nspec1d):
        re = datos.espectro.real[ii]
        im = datos.espectro.imag[ii]
        if autoph:
            spec1d, phase = ng.proc_autophase.autops(re+1j*im,
                                            "acme",

                                            return_phases=True,
                                            disp=False)
        else:
            spec1d = re 
            phase = [0]
        # spec1d_re = ng.proc_bl.cbf(spec1d.real)
        if absolute:
            spec1d_re = np.abs(spec1d)
        else:
            spec1d_re = spec1d.real
            # # Calculate baseline as a straight line between the mean of the first 100 and last 100 points
            # baseline_points = 50  # Number of points to use for baseline calculation
            # mean_start = np.mean(spec1d.real[:baseline_points])
            # mean_end = np.mean(spec1d.real[-baseline_points:])
            # mean_ppm_start = np.mean(ppmAxis[:baseline_points])
            # mean_ppm_end = np.mean(ppmAxis[-baseline_points:])
            # baseline = ppmAxis * (mean_end - mean_start) / (mean_ppm_end - mean_ppm_start) + mean_start
            # spec1d_re = spec1d.real-baseline  # remove DC offset
        spec_vs_t[ii, :] = spec1d_re
        phases.append(phase[0])
        ax_spec = axs_spec[jj]
        ax_spec.plot(ppmAxis, spec1d_re,
                    color=colors[jj], alpha=0.1)# 0.2+(ii/datos.espectro.size[0])*0.7)
        # find the ppm of the maximum in the range < ppm_treshold ppm
        re_in_ROI = spec1d_re # re in Region Of Interest
        ppmAxis_in_ROI = ppmAxis
        max_index = np.argmax(re_in_ROI)
        ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]
        ppm_of_max.append(ppm_of_max_in_equilibrium)
    ppm_of_max = np.array(ppm_of_max)
    spec_vs_t_list.append(spec_vs_t)


    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvline(r1, color='k', linestyle='--')
    ax_spec.axvline(r2, color='k', linestyle='--')
    ax_spec.set_xlabel('chemical shift [ppm]')
    ax_spec.set_ylabel('Intensity [a.u.]')
    ax_spec.axhline(0, color='k')


    ax.plot(tau, ppm_of_max, 'o-', color=colors[0])
    axph.plot(tau, phases, 'o-', color=colors[0])

    signal = datos.Integrar(ppmRange=ppmRange)
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


ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')
ax_int.set_xlim(0,90)

ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta_{max}$ [ppm]')
ax.set_xlim(0, 90)

axph.set_xlabel('Time [h]')
axph.set_ylabel('Phase [deg]')
axph.set_xlim(0, 90)

#%%
specs_vs_t_array = np.concatenate([arr for arr in spec_vs_t_list], axis=0)
fig, ax = plt.subplots(num=17856522, figsize=(3, 6))
ax.imshow(specs_vs_t_array, aspect='auto')
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$^7$LiChemical shift [ppm]')

# # guardo data:
# if save:
#     filename = f'{savepath}/{muestra}_T1.png'
#     fig.savefig(filename)   # save the figure to file

#     Signals = np.array(Signals).T
#     tau = tau.reshape(tau.size, 1)
#     T1data = np.hstack((tau, Signals))
#     header = "tau [s]\t"
#     for ppmRange in ppmRanges:
#         header += f"{ppmRange} ppm\t"
#     np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

#     data = np.array([ppmAxis, re, im]).T
#     np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)

# %%
