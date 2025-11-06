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




####### metal
# expns, sample = [3, 6,27], 'Li on Cu mesh'
# vdlists =["/lists/vd/sm2974.T1-32"]*3
# xlabels = [r'$\mu$w OFF - RT', r'$\mu$w OFF - LT', r'$\mu$w ON - LT']
# # expns = np.arange(31, 34)
# # expns = np.concatenate([np.arange(1, 22), np.arange(31, 34)])  # directorio de datos
# plotRange = [400,200]
# ppmIntegrationWidth = 70  # ancho de la ventana de integracion
# ppmIntegration = [320, 225]
# Integration_around_max = [True, False, True]

#####

expns, sample = [101, 103, 114], 'Li on Cu mesh'
vdlists =["/lists/vd/sm2974.T1-32","/lists/vd/sm2974.T1-16","/lists/vd/sm2974.T1-16"]
xlabels = [r'$\mu$w OFF - RT', r'$\mu$w OFF - LT', r'$\mu$w ON - LT']
# expns = np.arange(31, 34)
# expns = np.concatenate([np.arange(1, 22), np.arange(31, 34)])  # directorio de datos
plotRange = [600,200]
ppmIntegrationWidth = 70  # ancho de la ventana de integracion
ppmIntegration = [500, 400]
Integration_around_max = [True, True, False]
####

# for the initial fit:
center_gues = 280  # ppm
sigma_guess = 1  # ppm
gamma_guess = 1  # ppm
height_guess = 3.4e8  # a.u.


absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-10-27_InSitu/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-10_InSitu/"
muestra = ""


# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-06-26_InSitu/"
# expns = [37]
# plotRange = [500,0]
# muestra = "eLi-uwON"

# directorio de guradado

#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================

T1_list = []
T1_err_list = []
ppm_of_max_list = []  # lista de ppm del maximo de cada espectro
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
    datos = DatosProcesadosT1(f'{path}/{expn}/', vdlist_path=vdlists[jj])    
    datos.espectro.ppmSelect(plotRange)
    #+++++++++++++++++++++++++++++++
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    re = datos.espectro.real[-1]
    im = datos.espectro.imag[-1]

    if Integration_around_max[jj]:
        # Find the position of the maximum in re
        max_index = np.argmax(re)
        ppm_of_maximum = ppmAxis[max_index]
        print(f"Maximum value in re at: {ppm_of_maximum:.2f} ppm (index {max_index})")


        ppmRange = [ppm_of_maximum - ppmIntegrationWidth/2,
                    ppm_of_maximum + ppmIntegrationWidth/2]
    else:
        ppmRange = [min(ppmIntegration), max(ppmIntegration)]
    # region of interest
    # spec1d, phase = ng.proc_autophase.autops(re+1j*im,
    #                                         "acme",
    #                                         return_phases=True,
    #                                         disp=False)
    # re = spec1d.real
    # im = spec1d.imag
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    fig1d, ax1d = plt.subplots()
    title = f"expn: {expn}\n"
    ax1d.set_title(title)
    ax1d.axvspan(r1, r2, alpha=0.2, color='b')
    ax1d.axhline(0, color='gray')
    ax1d.plot(ppmAxis, re/np.max(re), 'k', lw=2)
    ax1d.text(r1-np.abs(0.1*r1), 0.8, "Region de integracion\n(T1)", color='b')
    ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax1d.set_xlabel("NMR Shift [ppm]")

    tau, signal = datos.get_T1data(ppmRange)
    tau_fit, signal_fit, residuals = datos.T1fit()

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(title)
    # -------------
    axs[0, 0].plot(tau, signal, 'ko')
    axs[0, 0].plot(tau_fit, signal_fit, 'b-')
    text = f"$T_1 =$ {datos.T1params[1]:.1e} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text, multialignment="left")
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    # -------------
    axs[1, 0].plot(tau, residuals, 'ko')
    axs[1, 0].axhline(0, color='k', linestyle='--')
    axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
    # -------------
    axs[0, 1].plot(tau, signal, 'ko')
    axs[0, 1].plot(tau_fit, signal_fit, 'b-')
    axs[0, 1].set(xlabel=r'$\tau$ [ms]')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_yticklabels([])
    # -------------
    axs[1, 1].plot(tau, residuals, 'ko')
    axs[1, 1].axhline(0, color='k', linestyle='--')
    axs[1, 1].set(xlabel=r'$\tau$ [ms]')
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_yticklabels([])
    

    T1_list.append(datos.T1params[1])  # guardo T1
    # T1_err_list.append(datos.T1params_err[1])
    T1_err_list.append(datos.T1stderr[1])


T1_list = np.array(T1_list)
ppm_of_max_list = np.array(ppm_of_max_list)
#%%
# Plot: T1 vs time
fig, ax = plt.subplots(figsize=(5, 4))
ax.set_title(sample)
bars = ax.bar(np.arange(T1_list.size), T1_list, yerr=T1_err_list, label='Experimental $T_1$', capsize=5)
# ax.set_xlabel("nexp")
ax.set_ylabel(r"$T_1$ [ms]")
ax.grid(axis='y')
ax.set_xticks(np.arange(T1_list.size))
ax.set_xticklabels(xlabels)
# add value labels above the bars
# if T1_list.size > 0:
#     pad = 0.02 * np.max(T1_list)  # vertical padding above bars
#     for bar in bars:
#         h = bar.get_height()
#         ax.text(bar.get_x() + bar.get_width() / 4,
#                 h + pad, f"{h:.0f}", ha='center', va='bottom')
# ax.set_ylim(0, 250)
#%%=====================================================================
## T from T1
# ## mw OFF calibration
# a = 43.7 *1000 # msK
# b = -3.2 * 1e-3# ms
# T1_s = T1_list * 1000
# T = a / (T1_s - b)
# print(T-T[0])

# ### mw ON calibration
# a = 40.9 # msK
# b = 10.6 * 1e-3 # ms
# T1_s = T1_list * 1000
# T = a / (T1_s - b)
# print(T-T[0])
# #%% Korringa

T_K = np.linspace(230, 305, 1024)
T1_K = 43.4 / T_K
plt.figure(figsize=(5, 4))
plt.plot(T_K, T1_K, 'k--', label='Korringa relation')
plt.xlabel('T [K]')
plt.ylabel(r'$T_1$ [s]')

T_from_T1 = 43.4 / (T1_list/1000)
T_err = np.abs(T_from_T1*(T1_err_list/T1_list))
plt.scatter(T_from_T1, T1_list/1000, color='r', label='Experimental data ')


#%%
plt.figure(figsize=(5, 4))
plt.title(sample)
# bars = plt.bar(np.arange(T1_list.size), T_from_T1, yerr=T_err, capsize=5)
# plt.ylabel(r"$T$ [K]")
bars = plt.bar(np.arange(T1_list.size), T_from_T1-T_from_T1[0], yerr=T_err, capsize=5)
plt.ylabel(r"$\Delta T$ [K]")
plt.grid(axis='y')

# add value labels above the bars
# if T_from_T1.size > 0:
#     pad = 0.01 * np.max(T_from_T1)  # vertical padding above bars
#     for bar, val in zip(bars, T_from_T1):
#         h = bar.get_height() 
#         plt.text(bar.get_x() + bar.get_width() / 2, h + pad, f"{val-273.15:.0f} "+r"$^{\circ}$C", ha='center', va='bottom')
# plt.ylim([0, 320])
plt.tight_layout()
plt.xticks(np.arange(T1_list.size), xlabels)
# %%
