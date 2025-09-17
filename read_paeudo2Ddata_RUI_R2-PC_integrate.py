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


# directorio de datos
expns = np.arange(30,234)
expns = np.arange(10, 201)

absolute= False
autoph = False
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-08-10_ccATMC_Rui-R1_LFP-Cu_PC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-08_R2/"
muestra = "7Li_cellR2-PCprotocol"

save = False
plotRange = [350, 150]
# rango de integracion
ppmRange = [300,200]

# #### diamagnetic
# plotRange = [100, -100]
# # rango de integracion
# ppmRange = [50,-50]


colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))

signals = np.zeros(expns.size)
tau = np.zeros(expns.size)
ppm_of_max = np.zeros(expns.size)
ppm_mean = np.zeros(expns.size)
phases = np.zeros(expns.size)
for jj, expn in enumerate(expns):
    print(f"Processing expn {expn}, {jj+1} of {expns.size}")
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    spec_time = datos.acqus.dic['DATE_START']
    # spec_time = datos.acqus.dic['DATE_START'] - Nominal_duration_of_experiment # to take into account the time of autotune
    if jj==0:
        t_0 = spec_time# s
        spec_vs_t = np.zeros([spec.size, expns.size])
    tau[jj] = (spec_time - t_0)/3600  # convert to hours

    # grafico todos los espectros juntos
    re = datos.espectro.real
    im = datos.espectro.imag
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
    spec_vs_t[:, jj] = spec1d_re
    phases[jj] = phase[0]
    ax_spec.plot(ppmAxis, spec1d_re)# 0.2+(ii/datos.espectro.size[0])*0.7)
    # find the ppm of the maximum in the range < ppm_treshold ppm
    re_in_ROI = spec1d_re # re in Region Of Interest
    ppmAxis_in_ROI = ppmAxis
    max_index = np.argmax(re_in_ROI)
    ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]
    ppm_of_max[jj] = ppm_of_max_in_equilibrium

    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvline(r1, color='k', linestyle='--')
    ax_spec.axvline(r2, color='k', linestyle='--')
    ax_spec.set_xlabel('chemical shift [ppm]')
    ax_spec.set_ylabel('Intensity [a.u.]')
    ax_spec.axhline(0, color='k')



    signal = datos.Integrar(ppmRange=ppmRange)
    if jj == 0:
        initial_signal = signal
    signals[jj] = signal  # normalizo la señal al primer espectro

    ppm_mean[jj] = datos.Mean_ppm(ppmRange=ppmRange)

#%% -------------
fig,ax = plt.subplots(num=1785731, figsize=(8, 3))
figph,axph = plt.subplots(num=178573111, figsize=(8, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(8, 3))


tau_real, tau_continuo = tau, tau

ax_int.plot(tau_continuo, signals/signals.max(), 'o')
ax_int.yaxis.set_major_locator(MultipleLocator(0.1))  # línea cada 0.1 en Y
ax_int.grid(axis='y')  # solo grilla en Y
ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')
#ax_int.set_xlim([-5, 70])

ax.plot(tau_continuo, ppm_of_max, 'o-', label='ppm of max signal')
ax.plot(tau_continuo, ppm_mean, 'o-', label='mean ppm in ROI')
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta_{max}$ [ppm]')
ax.legend()
ax.set_ylim(235.1, 259.9)

axph.plot(tau_continuo, phases, 'o-')
axph.set_xlabel('Time [h]')
axph.set_ylabel('Phase [deg]')
# axph.set_xlim(0, 90)


#%%
# Asegurarse de que sean 2D arrays de coordenadas
ppm_mesh, tau_mesh = np.meshgrid(ppmAxis, tau)

fig, ax = plt.subplots(num=17856522, figsize=(3, 6))
pcm = ax.pcolormesh(ppm_mesh, tau_mesh, spec_vs_t.T, shading='auto')  # shading='auto' ajusta la interpolación
ax.set_ylabel('Time [h]')
ax.set_xlabel(r'$^7$Li Chemical shift [ppm]')
fig.colorbar(pcm, ax=ax, label='Intensity')  # Opcional: barra de color
# ax.set_xlim(ppmRange)
ax.set_xlim([260, 220])
ax.axvline(ppm_of_max[0], color='k', linestyle='--')

# %%
