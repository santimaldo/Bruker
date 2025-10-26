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
expns = [8, 21]
samples = ["19F_LP57-neat", "19F_LP57-Mn(TFSI)2-8mM"]
plotRange = [100, -200]
ppmRange = [-63, -75] # ventana de integracion



# expns = [12, 24]
# samples = ["7Li_LP57-neat", "7Li_LP57-Mn(TFSI)2-8mM"]
# plotRange = [100, -100]
# ppmIntegrationWidth = 20  # ancho de la ventana de integracion

absolute= False
autoph = False
save = False

path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-10-13_DRinsitu_LP57-MnTFSI/"
# directorio de guradado
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\TMions\Analysis\2025-10_in-situ_relaxation/"
muestra = "eLi-KBr"

T1_short_list = []
T1_long_list = []
A_short_list = []
A_long_list = []
#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
# Obtener la paleta Set1
cmap = plt.get_cmap("Set1")
colors = [cmap(i) for i in range(cmap.N)]  # 9 colores
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, [expn, sample] in enumerate(zip(expns, samples)):
    color=colors[jj]
    #=====================================================================
    # Ajuste de espectro 2D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesadosT1(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)
    #+++++++++++++++++++++++++++++++
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    re = datos.espectro.real[-1]
    im = datos.espectro.imag[-1]

    # Find the position of the maximum in re
    # max_index = np.argmax(re)
    # ppm_of_maximum = ppmAxis[max_index]
    # print(f"Maximum value in re at: {ppm_of_maximum:.2f} ppm (index {max_index})")

    # ppmRange = [ppm_of_maximum - ppmIntegrationWidth/2,
    #             ppm_of_maximum + ppmIntegrationWidth/2]  # region of interest
    # spec1d, phase = ng.proc_autophase.autops(re+1j*im,
    #                                         "acme",

    #                                         return_phases=True,
    #                                         disp=False)
    # re = spec1d.real
    # im = spec1d.imag
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    fig1d, ax1d = plt.subplots()
    title = f"sample: {sample}"
    ax1d.set_title(title)
    ax1d.axvspan(r1, r2, alpha=0.2, color=color)
    ax1d.axhline(0, color='gray')
    ax1d.plot(ppmAxis, re/np.max(re), 'k', lw=2)
    ax1d.text(r1-np.abs(0.1*r1), 0.8, "Region de integracion\n(T1)", color=color)
    ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax1d.set_xlabel("NMR Shift [ppm]")

    tau, signal = datos.get_T1data(ppmRange)
    tau_fit, signal_fit, residuals = datos.T1fit(model='bi')
    T1_1, T1_2 = datos.T1params[1], datos.T1params[3]
    T1_short = min([T1_1, T1_2])
    T1_long = max([T1_1, T1_2])
    T1_short_list.append(T1_short / 1000)  # convertir de ms a s
    T1_long_list.append(T1_long / 1000)  # convertir de ms a s
    
    A1, A2, y0 = datos.T1params[0], datos.T1params[2], datos.T1params[4]
    c1 = A1 / (A1 + A2)
    c2 = A2 / (A1 + A2)
    if T1_short == T1_1:
        A_short = c1
        A_long = c2
    else:
        A_short = c2
        A_long = c1
    A_short_list.append(A_short)
    A_long_list.append(A_long)    

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(title)
    # -------------
    axs[0, 0].plot(tau, signal, 'ko')
    axs[0, 0].plot(tau_fit, signal_fit, '-', color=color)
    text = f"$T_1 =$ {T1_short:.0f} ms  ({A_short*100:.0f} %) \n $T_1 =$ {T1_long:.0f} ms  ({A_long*100:.0f} %) "
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
                multialignment="left")
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    # -------------
    axs[1, 0].plot(tau, residuals, 'ko')
    axs[1, 0].axhline(0, color='k', linestyle='--')
    axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
    # -------------
    axs[0, 1].plot(tau, signal, 'ko')
    axs[0, 1].plot(tau_fit, signal_fit, '-', color=color)
    axs[0, 1].set(xlabel=r'$\tau$ [ms]')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_yticklabels([])
    # -------------
    axs[1, 1].plot(tau, residuals, 'ko')
    axs[1, 1].axhline(0, color='k', linestyle='--')
    axs[1, 1].set(xlabel=r'$\tau$ [ms]')
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_yticklabels([])







#%%
print("\nResumen de T1 y temperaturas estimadas:")
print("Expn\tT1 [s]\tT_off [K]\tT_on [K]")

for expn, T1 in zip(expns, T1_list):
    try:
        T_off = 43.7 / T1  - 3.2  # T1 en s
    except ZeroDivisionError:
        T_off = np.nan
    try:
        T_on = 40.9 / T1 + 10.6 # T1 en s 
    except ZeroDivisionError:
        T_on = np.nan
    print(f"{expn}\t{T1:.3f}\t{T_off:.1f}\t\t{T_on:.1f}")

# %%
