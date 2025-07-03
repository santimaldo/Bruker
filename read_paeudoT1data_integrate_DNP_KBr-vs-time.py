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
expns = np.concatenate([np.arange(1, 60), np.arange(61, 100)])  # directorio de datos

plotRange = [60, 10]
ppmIntegrationWidth = 40  # ancho de la ventana de integracion

# for the initial fit:
center_gues = 36.2  # ppm
sigma_guess = 1  # ppm
gamma_guess = 1  # ppm
height_guess = 3.4e8  # a.u.


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

times_list = []  # lista de tiempos de inicio de cada experimento
T1_list = []
ppm_of_max_list = []  # lista de ppm del maximo de cada espectro
colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
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
    max_index = np.argmax(re)
    ppm_of_maximum = ppmAxis[max_index]
    print(f"Maximum value in re at: {ppm_of_maximum:.2f} ppm (index {max_index})")

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


    vfit = vf.VoigtFit(ppmAxis, re, Npicos=1, center=center_gues, sigma=sigma_guess,
                      gamma=gamma_guess, height=height_guess)
    fig_vfit = vfit.plot_ajuste()
    ax_vfit = fig_vfit.gca()
    ax_vfit.set_xlim(ppmRange)
    ax_vfit.set_title(title)
    # center_guess = [vfit.params["m1_center"].value]
    # sigma_guess = [vfit.params["m1_sigma"].value]
    # gamma_guess = [vfit.params["m1_gamma"].value]
    # height_guess = [vfit.params["m1_amplitude"].value]
    ppm_of_max_list.append(vfit.params["m1_center"].value)  # guardo el ppm del maximo
 
    tau, signal = datos.get_T1data(ppmRange)
    tau_fit, signal_fit, residuals = datos.T1fit()


    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(title)
    # -------------
    axs[0, 0].plot(tau, signal, 'ko')
    axs[0, 0].plot(tau_fit, signal_fit, 'b-')
    text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
                multialignment="left")
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

# --------------------------------------
# Thurber empirical model for T1 vs T
# --------------------------------------
def thurber_model(T):
    """
    Empirical model by Thurber for T1 relaxation time as a function of temperature (K).
    """
    return 0.0145 + 5330 * T**-2 + 1.42e7 * T**-4 + 2.48e9 * T**-6

# Plot: T1 vs time
plt.figure(figsize=(6, 4))
plt.plot(time_minutes, T1_list, 'o-', label='Experimental $T_1$')
plt.xlabel("Time [min]")
plt.ylabel(r"$T_1$ [ms]")
plt.title(r"$T_1$ as a function of time")
plt.grid(True)
plt.tight_layout()

# Inverse interpolation: T1 → T
T = np.linspace(290, 350, 1000)
T1_model = thurber_model(T)
T1_to_T = interp1d(T1_model, T, kind='linear', fill_value='extrapolate')

# Temperature estimated from T1
T_from_T1_C = T1_to_T(T1_list / 1000) - 273.15  # in °C

plt.figure(figsize=(6, 4))
plt.plot(time_minutes, T_from_T1_C, 'o-', label='Estimated T from $T_1$')
plt.xlabel("Time [min]")
plt.ylabel("Estimated temperature [°C]")
plt.title("Temperature estimated from $T_1$")
plt.grid(True)
plt.tight_layout()

# Plot: T1 vs Temperature (model vs experimental)
fig, ax = plt.subplots(figsize=(8, 4))
T_plot = np.linspace(20, 296, 1000)  # Temperature range fixed as requested
T1_plot = thurber_model(T_plot)

ax.plot(T_plot, T1_plot, label='Thurber model')
ax.plot(T1_to_T(T1_list / 1000), T1_list / 1000, 'o-', label='Experimental data')
ax.set_xlabel("Temperature [K]")
ax.set_ylabel(r"$T_1$ [s]")
ax.set_yscale('log')
ax.set_title(r"$T_1$ vs Temperature")
ax.legend()
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
fig.tight_layout()

# --- Chemical shift data cleaning ---
print("WARNING: Removing points with chemical shift > 37 ppm")
valid_idx = ppm_of_max_list < 37
t_valid = time_minutes[valid_idx]
ppm_valid = np.array(ppm_of_max_list)[valid_idx]

plt.figure(figsize=(6, 4))
plt.plot(t_valid, ppm_valid, 'o-', label='KBr chemical shift')
plt.xlabel("Time [min]")
plt.ylabel("Chemical shift [ppm]")
plt.title("KBr chemical shift over time")
plt.grid(True)
plt.tight_layout()

# --- Chemical shift to temperature calibration ---
T_chemshift_slope = -0.025  # ppm/K
T0 = T1_to_T(T1_list[0] / 1000)
ppm0 = ppm_valid[0]
offset_ppm = ppm0 - T_chemshift_slope * T0

def chemshift_to_T(ppm):
    """
    Convert chemical shift (ppm) to temperature (K).
    """
    return (ppm - offset_ppm) / T_chemshift_slope

# Estimated temperature from chemical shift
T_from_shift_C = chemshift_to_T(ppm_valid) - 273.15

plt.figure(figsize=(6, 4))
plt.plot(t_valid, T_from_shift_C, 'o-', label='Estimated T from chemical shift')
plt.xlabel("Time [min]")
plt.ylabel("Estimated temperature [°C]")
plt.title("Temperature estimated from chemical shift")
plt.grid(True)
plt.tight_layout()

# --- Combined plot: temperature estimates from both methods ---
plt.figure(figsize=(7, 4))
plt.plot(time_minutes, T_from_T1_C, 'o-', label='From $T_1$')
plt.plot(t_valid, T_from_shift_C, 'o-', label='From chemical shift')
plt.xlabel("Time [min]")
plt.ylabel("Estimated temperature [°C]")
plt.title("Comparison of temperature estimation methods")
plt.legend()
plt.grid(True)
plt.tight_layout()




#%%
    #++++++++++++++++++++++++++++    
    # time_counter = datos.acqus.dic['DATE_START']    
    # tau = np.linspace(t_ini, t_end, Nspec1d) - t_0 # horas
    # tau = tau/3600  # convert to hours
    # # grafico todos los espectros juntos
    # ppm_of_max = []
    # phases = []
    # for ii in range(Nspec1d):
    #     re = datos.espectro.real[ii]
    #     im = datos.espectro.imag[ii]
    #     if autoph:
    #         spec1d, phase = ng.proc_autophase.autops(re+1j*im,
    #                                         "acme",

    #                                         return_phases=True,
    #                                         disp=False)
    #     else:
    #         spec1d = re 
    #         phase = [0]
    #     # spec1d_re = ng.proc_bl.cbf(spec1d.real)
    #     if absolute:
    #         spec1d_re = np.abs(spec1d)
    #     else:
    #         spec1d_re = spec1d.real
    #         # # Calculate baseline as a straight line between the mean of the first 100 and last 100 points
    #         # baseline_points = 50  # Number of points to use for baseline calculation
    #         # mean_start = np.mean(spec1d.real[:baseline_points])
    #         # mean_end = np.mean(spec1d.real[-baseline_points:])
    #         # mean_ppm_start = np.mean(ppmAxis[:baseline_points])
    #         # mean_ppm_end = np.mean(ppmAxis[-baseline_points:])
    #         # baseline = ppmAxis * (mean_end - mean_start) / (mean_ppm_end - mean_ppm_start) + mean_start
    #         # spec1d_re = spec1d.real-baseline  # remove DC offset
    #     spec_vs_t[ii, :] = spec1d_re
    #     phases.append(phase[0])
    #     ax_spec = axs_spec[jj]
    #     ax_spec.plot(ppmAxis, spec1d_re,
    #                 color=colors[jj], alpha=0.1)# 0.2+(ii/datos.espectro.size[0])*0.7)
    #     # find the ppm of the maximum in the range < ppm_treshold ppm
    #     re_in_ROI = spec1d_re # re in Region Of Interest
    #     ppmAxis_in_ROI = ppmAxis
    #     max_index = np.argmax(re_in_ROI)
    #     ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]
    #     ppm_of_max.append(ppm_of_max_in_equilibrium)
    # ppm_of_max = np.array(ppm_of_max)
    # spec_vs_t_list.append(spec_vs_t)


#     ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
#     r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
#     ax_spec.axvline(r1, color='k', linestyle='--')
#     ax_spec.axvline(r2, color='k', linestyle='--')
#     ax_spec.set_xlabel('chemical shift [ppm]')
#     ax_spec.set_ylabel('Intensity [a.u.]')
#     ax_spec.axhline(0, color='k')


#     ax.plot(tau, ppm_of_max, 'o-', color=colors[0])
#     axph.plot(tau, phases, 'o-', color=colors[0])

#     signal = datos.Integrar(ppmRange=ppmRange)
#     if jj == 0:
#         initial_signals = signal[0]
#     signal = signal/initial_signals  # normalizo la señal al primer espectro

#     # -------------
#     ax_int.plot(tau, signal[:tau.size], 'o', color=colors[0])
#     if expn==5:
#         # firts plating
#         fpd = 0.8 # first plating duration in hours
#         condition = tau<tau[0]+fpd
#         tau_plating = tau[condition]

#         signal_plating = signal[condition]
#         ax_int.plot(tau_plating, signal_plating, 'o', color=colors[jj])
#         ppm_of_max_plating = ppm_of_max[condition]
#         ax.plot(tau_plating, ppm_of_max_plating, 'o-', color=colors[jj])


# ax_int.set_xlabel('Time [h]')
# ax_int.set_ylabel('Normalized area')
# ax_int.set_xlim(0,90)

# ax.set_xlabel('Time [h]')
# ax.set_ylabel(r'$\delta_{max}$ [ppm]')
# ax.set_xlim(0, 90)

# axph.set_xlabel('Time [h]')
# axph.set_ylabel('Phase [deg]')
# axph.set_xlim(0, 90)

# #%%
# specs_vs_t_array = np.concatenate([arr for arr in spec_vs_t_list], axis=0)
# fig, ax = plt.subplots(num=17856522, figsize=(3, 6))
# ax.imshow(specs_vs_t_array, aspect='auto')
# ax.set_xlabel('Time [h]')
# ax.set_ylabel(r'$^7$LiChemical shift [ppm]')

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
