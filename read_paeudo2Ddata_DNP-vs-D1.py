# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 2025

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate
import re
from redor import redor_curve

def label_curve(ax, x, y, label, idx, offset=(0, 0), **kwargs):
    """Anota una curva en el punto `idx`, sin rotación."""
    x0 = x[idx] + offset[0]
    y0 = y[idx] + offset[1]

    ax.text(x0, y0, label,
            rotation=0,
            ha='center', va='center',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.5),
            **kwargs)

############################################################

# # Directorio de datos
# expn = 44
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-06-21_PEO-solid-electrolyte/"
# # Directorio de guardado
# savepath = r"C:/"
# muestra = ""
# save = False
# plot_individual_pairs = False  # Activar/desactivar gráficos por par
# plotRange = [-40, -100]

# Directorio de datos
expns = np.arange(10,27)
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-13_3.2mm_Rui-dendrites_vsD1/"
# Directorio de guardado
savepath = r"C:/"
muestra = ""
save = False
plot_individual_pairs = False  # Activar/desactivar gráficos por par
plotRange = [50, -10]




# Rango de integración
ppmRange = [40, 7]


#=====================================================================
#=====================================================================
Signals = np.zeros(len(expns))
Signals_OFF = np.zeros(len(expns))
D1s = np.zeros(len(expns))
for jj, expn in enumerate(expns):
    path_ON = f"{path}/{expn}/"
    datos = DatosProcesados(path_ON, read_pp=False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    D1s[jj] = datos.acqus.D1

    path_OFF = f"{path}/{expn+100}/"
    datos_OFF = DatosProcesados(path_OFF, read_pp=False)
    datos_OFF.espectro.ppmSelect(plotRange)
    spec_OFF = datos_OFF.espectro.real
    ppmAxis_OFF = datos_OFF.espectro.ppmAxis

    # Integrar cada espectro y analizar
    fig_spec, ax_spec = plt.subplots(num=expn, figsize=(8,5))
    colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
    ax_spec.plot(ppmAxis, spec.T, "r", label='ON')
    ax_spec.plot(ppmAxis_OFF, spec_OFF.T, "b", label='OFF')
    ax_spec.set_xlabel("ppm")
    
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]
    ax_spec.axvspan(r1, r2, alpha=0.15, color='grey')
    signal = datos.Integrar(ppmRange=ppmRange)
    signal_OFF = datos_OFF.Integrar(ppmRange=ppmRange)
    Signals[jj] = signal
    Signals_OFF[jj] = signal_OFF

#%%
fig, ax = plt.subplots(figsize=(8,5))
ax.plot(D1s, Signals, 'ro-', label='ON')
ax.plot(D1s, Signals_OFF, 'bo-', label='OFF')
ax.set_xlabel("Recycling Delay (s)")
ax.set_ylabel("Signal (a.u.)")
# ax.set_xscale("log")
ax.legend()

fig_ratio,  ax_ratio = plt.subplots(figsize=(8,5))
ax_ratio.plot(D1s, Signals/Signals_OFF, 'ko-')
ax_ratio.set_xlabel("Recycling Delay (s)")
ax_ratio.set_ylabel("DNP Enhancement")
# ax_ratio.set_xscale("log")
plt.show()
# %%
