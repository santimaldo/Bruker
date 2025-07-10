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

#%%=========================================================
def nmr_time_to_cycle_time(t):
    """
    Convierte el tiempo NMR a:
    - tiempo recortado (dentro del intervalo válido más cercano),
    - tiempo continuo (pegando intervalos válidos sin huecos, manteniendo el que contiene t=0 sin corrimiento).

    Parámetros:
        t (float or array-like): Tiempo(s) en escala absoluta.

    Returns:
        t_trimmed (np.array): Tiempos recortados (valores válidos dentro de los ciclos definidos).
        t_continuous (np.array): Tiempos ajustados a un eje continuo sin huecos.
    """
    time_limits = [[-0.6354857804192499, -0.11170973809539039],
                   [-0.074153277778116, 2.9192848376792],
                   [3.1725133888885506, 3.217035165541614],
                   [4.62957266697523, 4.6814746656640756],
                   [5.835546000328407, 6.797969031571111]]

    # Asegurar que t sea array
    t = np.atleast_1d(t)

    # --- Paso 1: encontrar el intervalo que contiene t = 0 ---
    zero_index = None
    for i, (start, end) in enumerate(time_limits):
        if start <= 0 <= end:
            zero_index = i
            break
    if zero_index is None:
        raise ValueError("t = 0 no está contenido en ningún intervalo de time_limits")

    # --- Paso 2: calcular offsets para cada intervalo ---
    offsets = [0.0] * len(time_limits)
    
    # Hacia adelante (después de t = 0)
    for i in range(zero_index + 1, len(time_limits)):
        prev_end = time_limits[i - 1][1]
        curr_start = time_limits[i][0]
        offsets[i] = offsets[i - 1] + (prev_end - curr_start)
    
    # Hacia atrás (antes de t = 0)
    for i in reversed(range(0, zero_index)):
        next_start = time_limits[i + 1][0]
        curr_end = time_limits[i][1]
        offsets[i] = offsets[i + 1] - (curr_end - next_start)

    # --- Paso 3: convertir t ---
    t_trimmed = []
    t_continuous = []

    for val in t:
        for i, (start, end) in enumerate(time_limits):
            if start <= val <= end:
                # Está dentro de un intervalo válido
                t_trimmed.append(val)
                t_continuous.append(val + offsets[i])
                break
            elif val < start:
                # Está en un hueco anterior a este intervalo
                if i == 0:
                    t_trimmed.append(start)
                    t_continuous.append(start + offsets[i])
                else:
                    t_trimmed.append(time_limits[i - 1][1])
                    t_continuous.append(time_limits[i - 1][1] + offsets[i - 1])
                break
        else:
            # Está después del último intervalo
            t_trimmed.append(val)
            t_continuous.append(val + offsets[-1])
    return np.array(t_trimmed), np.array(t_continuous)
#%%=========================================================


# directorio de datos
expns = np.arange(11,310)

Nominal_duration_of_experiment = 430 # seconds (7 minutes 10 seconds)
absolute= False
autoph = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-07-08_in-situ_IMEC-overlaid/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\in-situ\2025-07_eLI-Li_overlaid/"
muestra = "7Li_eLi-Li-cell_LP40"

save = False
plotRange = [400, 100]
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
phases = np.zeros(expns.size)
for jj, expn in enumerate(expns):
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    # spec_time = datos.acqus.dic['DATE_START']
    spec_time = datos.acqus.dic['DATE_START'] - Nominal_duration_of_experiment # to take into account the time of autotune
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

#%% -------------
fig,ax = plt.subplots(num=1785731, figsize=(8, 3))
figph,axph = plt.subplots(num=178573111, figsize=(8, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(8, 3))
fig_int.suptitle(muestra)


tau_real, tau_continuo = nmr_time_to_cycle_time(tau)  # convert NMR time to cycle time

ax_int.plot(tau_continuo, signals/signals[0], 'o')
ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')

ax.plot(tau_continuo, ppm_of_max, 'o-')
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta_{max}$ [ppm]')
# ax.set_xlim(0, 90)

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
ax.set_xlim([280, 240])
ax.axvline(ppm_of_max[0], color='k', linestyle='--')

# %%
