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
from VoigtFit import VoigtFit
import pandas as pd


# directorio de datos
expns = np.arange(30,234)
expns = np.arange(223, 221, -1)

absolute= False
autoph = False
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-08-08_ccATMC_Rui-R1_LFP-Cu_CC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-08_R1/"
muestra = "7Li_cellR1-CCprotocol"

save = False
plotRange = [350, 150]
# rango de integracion
ppmRange = [300,200]

peaks = [245, 260]
range_of_peaks_to_save = [-50, -65]  # Rango de ppm para guardar los picos

# initial guesses for Voigt fit
m1_amplitude, m2_amplitude = [4479075, 34970261]
m1_center, m2_center =  [239.5, 245]
m1_sigma, m2_sigma = [4.17, 7.5]

vfit_results = pd.DataFrame(columns=['expn', 'time',
                                  'm1_amplitude', 'm1_amplitude_stderr',
                                  'm1_center', 'm1_center_stderr',
                                  'm1_sigma', 'm1_sigma_stderr',
                                  'm2_amplitude', 'm2_amplitude_stderr',
                                  'm2_center', 'm2_center_stderr',
                                  'm2_sigma', 'm2_sigma_stderr'
                                  ])

colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))
for jj, expn in enumerate(expns):
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    spec_time = datos.acqus.dic['DATE_START']

    # grafico todos los espectros juntos
    spec1d_re = datos.espectro.real
    if jj==0:
        t_0 = spec_time# s
        spec_vs_t = np.zeros([spec.size, expns.size])
    spec_vs_t[:, jj] = spec1d_re
    ax_spec.plot(ppmAxis, spec1d_re)# 0.2+(ii/datos.espectro.size[0])*0.7)
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvline(r1, color='k', linestyle='--')
    ax_spec.axvline(r2, color='k', linestyle='--')
    ax_spec.set_xlabel('chemical shift [ppm]')
    ax_spec.set_ylabel('Intensity [a.u.]')
    ax_spec.axhline(0, color='k')

    ### ahora si, el ajuste de Voigt
    ydata = spec1d_re
    xdata = ppmAxis
    vfit=VoigtFit(xdata, 
              ydata, 
              Npicos=2,
              ajustar=True,
              amplitude=[m1_amplitude, m2_amplitude],
              center=[m1_center, m2_center],
              sigma=[m1_sigma, m2_sigma],
              #fijar=['center', 'sigma', 'gamma'],
              )
    fig = vfit.plot_ajuste()
    fig.gca().set_title(f"expn. {expn}")

    # redefino para ser usado en la proxima vuelta
    m1_amplitude = vfit.params["m1_amplitude"].value
    m2_amplitude = vfit.params["m2_amplitude"].value
    m1_center = vfit.params["m1_center"].value
    m2_center = vfit.params["m2_center"].value
    m1_sigma = vfit.params["m1_sigma"].value
    m2_sigma = vfit.params["m2_sigma"].value

    df = pd.DataFrame({'expn': [expn],
                      'time': spec_time,
                      'm1_amplitude': [vfit.params["m1_amplitude"].value],
                      'm1_amplitude_stderr': [vfit.params["m1_amplitude"].stderr],
                      'm1_center': [vfit.params["m1_center"].value],
                      'm1_center_stderr': [vfit.params["m1_center"].stderr],
                      'm1_sigma': [vfit.params["m1_sigma"].value],
                      'm1_sigma_stderr': [vfit.params["m1_sigma"].stderr],
                      'm2_amplitude': [vfit.params["m2_amplitude"].value],
                      'm2_amplitude_stderr': [vfit.params["m2_amplitude"].stderr],
                      'm2_center': [vfit.params["m2_center"].value],
                      'm2_center_stderr': [vfit.params["m2_center"].stderr],
                      'm2_sigma': [vfit.params["m2_sigma"].value],
                      'm2_sigma_stderr': [vfit.params["m2_sigma"].stderr],
                      })
    vfit_results = pd.concat([vfit_results, df])

#%% -------------
fig,ax = plt.subplots(num=1785731, figsize=(8, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(8, 3))

time = vfit_results['time'] - vfit_results['time'].min()
time = time/3600  # time in hours

signals = vfit_results['m1_amplitude'] + vfit_results['m2_amplitude']
# signals_stderr = np.sqrt(vfit_results['m1_amplitude_stderr']**2 + vfit_results['m2_amplitude_stderr']**2)
# signals_err = 1.96*signals_stderr  # 95% confidence interval

ax_int.plot(time, vfit_results['m1_amplitude']/signals.max(), 'o', label='peak 1')
ax_int.plot(time, vfit_results['m2_amplitude']/signals.max(), 'o', label='peak 2')
ax_int.plot(time, signals/signals.max(), 'ko', label='total signal')
ax_int.yaxis.set_major_locator(MultipleLocator(0.1))  # línea cada 0.1 en Y
ax_int.grid(axis='y')  # solo grilla en Y
ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')
#ax_int.set_xlim([-5, 70])

ax.plot(time, vfit_results['m1_center'], 'o-', label='peak 1')
ax.plot(time, vfit_results['m2_center'], 'o-', label='peak 2')
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta_{max}$ [ppm]')
ax.legend()
# ax.set_xlim(0, 90)



#%%
# Asegurarse de que sean 2D arrays de coordenadas
ppm_mesh, tau_mesh = np.meshgrid(ppmAxis, time)

fig, ax = plt.subplots(num=17856522, figsize=(3, 6))
pcm = ax.pcolormesh(ppm_mesh, tau_mesh, spec_vs_t.T, shading='auto')  # shading='auto' ajusta la interpolación
ax.set_ylabel('Time [h]')
ax.set_xlabel(r'$^7$Li Chemical shift [ppm]')
fig.colorbar(pcm, ax=ax, label='Intensity')  # Opcional: barra de color
# ax.set_xlim(ppmRange)
ax.set_xlim([260, 220])
# ax.axvline(ppm_of_max[0], color='k', linestyle='--')

# %%
