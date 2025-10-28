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
# expns = np.concatenate([np.arange(1, 60), np.arange(61, 100)])  # directorio de datos
expns = np.arange(2, 16)
# expns = np.arange(31, 34)
# expns = np.concatenate([np.arange(1, 22), np.arange(31, 34)])  # directorio de datos
plotRange = [500,100]
ppmIntegrationWidth = 150  # ancho de la ventana de integracion


# expns = np.arange(40, 55)
# plotRange = [500, 300]
# ppmIntegrationWidth = 60  # ancho de la ventana de integracion


# for the initial fit:
center_guess = 290 # ppm
sigma_guess = 1  # ppm
gamma_guess = 1  # ppm
height_guess = 3.4e8  # a.u.


absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-10-27_InSitu-T1/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-10_InSitu\Li-CuMesh-T1-while-cooling-down/"
muestra = "Li on Cu Mesh"
vdlist_path = r"/lists\vd\sm2974.T1-16" # relativo


# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-06-26_InSitu/"
# expns = [37]
# plotRange = [500,0]
# muestra = "eLi-uwON"

# directorio de guradado

#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================

times_list = []  # lista de tiempos de inicio de cada experimento
T1_list = []
ppm_of_max_list = []  # lista de ppm del maximo de cada espectro
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
    datos = DatosProcesadosT1(f'{path}/{expn}/', vdlist_path=vdlist_path)    
    datos.espectro.ppmSelect(plotRange)
    #+++++++++++++++++++++++++++++++
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    re = datos.espectro.real[-1]
    im = datos.espectro.imag[-1]

    # Find the position of the maximum in re
    max_index = np.argmax(re)
    ppm_of_maximum = ppmAxis[max_index]
    #### en lugar del maximo, uso el centro de masa
    #ppm_of_maximum = np.trapezoid(re*ppmAxis, x=ppmAxis)/np.trapezoid(re, x=ppmAxis)
    
    ppmRange = [ppm_of_maximum - ppmIntegrationWidth/2,
                ppm_of_maximum + ppmIntegrationWidth/2]  # region of interest
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
    text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
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
    
    times_list.append(datos.acqus.dic['DATE_START'])  # guardo el tiempo de inicio del experimento
    T1_list.append(datos.T1params[1])  # guardo T1




times_list = np.array(times_list)
T1_list = np.array(T1_list)
ppm_of_max_list = np.array(ppm_of_max_list)
#%%
# Convert time to minutes, starting from zero
time_minutes = (np.array(times_list) - np.min(times_list)) / 60


# Plot: T2 vs time
plt.figure(figsize=(6, 4))
plt.plot(time_minutes, T1_list, 'o-', label='Experimental $T_1$')
plt.xlabel("Time [min]")
plt.ylabel(r"$T_1$ [ms]")
plt.title(r"$T_1$ as a function of time")
plt.grid(True)
plt.tight_layout()
plt.xlim([-2, 79])

# # Inverse interpolation: T1 → T
# T = np.linspace(290, 350, 1000)
# T1_model = thurber_model(T)
# T1_to_T = interp1d(T1_model, T, kind='linear', fill_value='extrapolate')

# # Temperature estimated from T1
# T_from_T1_C = T1_to_T(T1_list / 1000) - 273.15  # in °C

# plt.figure(figsize=(6, 4))
# plt.plot(time_minutes, T_from_T1_C, 'o-', label='Estimated T from $T_1$')
# plt.xlabel("Time [min]")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Temperature estimated from $T_1$")
# plt.grid(True)
# plt.tight_layout()

# # Plot: T1 vs Temperature (model vs experimental)
# fig, ax = plt.subplots(figsize=(8, 4))
# T_plot = np.linspace(20, 296, 1000)  # Temperature range fixed as requested
# T1_plot = thurber_model(T_plot)

# ax.plot(T_plot, T1_plot, label='Thurber model')
# ax.plot(T1_to_T(T1_list / 1000), T1_list / 1000, 'o-', label='Experimental data')
# ax.set_xlabel("Temperature [K]")
# ax.set_ylabel(r"$T_1$ [s]")
# ax.set_yscale('log')
# ax.set_title(r"$T_1$ vs Temperature")
# ax.legend()
# ax.grid(True, which='both', linestyle='--', linewidth=0.5)
# fig.tight_layout()

# # --- Chemical shift data cleaning ---
# print("WARNING: Removing points with chemical shift > 37 ppm")
# valid_idx = ppm_of_max_list < 37
# t_valid = time_minutes[valid_idx]
# ppm_valid = np.array(ppm_of_max_list)[valid_idx]

# plt.figure(figsize=(6, 4))
# plt.plot(t_valid, ppm_valid, 'o-', label='KBr chemical shift')
# plt.xlabel("Time [min]")
# plt.ylabel("Chemical shift [ppm]")
# plt.title("KBr chemical shift over time")
# plt.grid(True)
# plt.tight_layout()

# # --- Chemical shift to temperature calibration ---
# T_chemshift_slope = -0.025  # ppm/K
# T0 = T1_to_T(T1_list[0] / 1000)
# ppm0 = ppm_valid[0]
# offset_ppm = ppm0 - T_chemshift_slope * T0

# def chemshift_to_T(ppm):
#     """
#     Convert chemical shift (ppm) to temperature (K).
#     """
#     return (ppm - offset_ppm) / T_chemshift_slope

# # Estimated temperature from chemical shift
# T_from_shift_C = chemshift_to_T(ppm_valid) - 273.15

# plt.figure(figsize=(6, 4))
# plt.plot(t_valid, T_from_shift_C, 'o-', label='Estimated T from chemical shift')
# plt.xlabel("Time [min]")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Temperature estimated from chemical shift")
# plt.grid(True)
# plt.tight_layout()

# # --- Combined plot: temperature estimates from both methods ---
# plt.figure(figsize=(7, 4))
# plt.plot(time_minutes, T_from_T1_C, 'o-', label='From $T_1$')
# plt.plot(t_valid, T_from_shift_C, 'o-', label='From chemical shift')
# plt.xlabel("Time [min]")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Comparison of temperature estimation methods")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()



# %%
