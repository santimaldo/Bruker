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

# Directorio de datos
expn = 44
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-06-21_PEO-solid-electrolyte/"
# Directorio de guardado
savepath = r"C:/"
muestra = ""
save = False
plot_individual_pairs = False  # Activar/desactivar gráficos por par
plotRange = [-40, -100]

# Rango de integración
ppmRanges = [
    [-57, -62],
    [-80, -81]
]

#=====================================================================
# Lectura del experimento 2D
#=====================================================================

path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(path_2D, read_pp=False)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

# Graficar los espectros 1D
fig_spec, ax_spec = plt.subplots(num=382910)
for kk in range(spec.shape[0]):
    ax_spec.plot(ppmAxis, spec[kk, :])
ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax_spec.axhline(0, color='k')

# Integrar cada espectro y analizar
colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
for ii, ppmRange in enumerate(ppmRanges):
    Signals = np.array([])
    color = colors[ii % len(colors)]
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]
    ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
    signal = datos.Integrar(ppmRange=ppmRange)
    Signals = np.append(Signals, signal)

    #=====================================================================
    # Calculo y grafico de (S - S0)/S0
    #=====================================================================

    # Separar S y S0
    S = Signals[0::2]
    S0 = Signals[1::2]
    N = np.arange(1, len(S) + 1) # number of rotor cycles

    # Calcular la razón (S - S0)/S0
    delta_S = (S0 - S) / S0
    spin_speed = 14000  # spinning speed in Hz
    recopl_time = N / spin_speed  # recoupling time in seconds

    # Graficar el resultado
    fig_redor, ax_redor = plt.subplots()
    ax_redor.plot(delta_S, 'o-')
    ax_redor.set_xlabel("Índice (simula tiempo)")
    ax_redor.set_ylabel("(S - S0) / S0")
    ax_redor.set_title(f"Redor in range: {ppmRange[0]} to {ppmRange[1]} ppm")
    ax_redor.grid(True)
    ax_redor.set_ylim(-0.1, 1.1)

    # Graficar S y S0 en bruto
    fig_t2, ax_t2 = plt.subplots()
    ax_t2.plot(S, 'o-', label='S')
    ax_t2.plot(S0, 'o-', label='S0')
    ax_t2.set_xlabel("Índice (simula tiempo)")
    ax_t2.set_ylabel("Señal integrada")
    ax_t2.set_title(f"ppm range: {ppmRange[0]} to {ppmRange[1]}")
    ax_t2.grid(True)
    ax_t2.legend()


    # Graficar S y S0 en bruto
    fig_redor2, ax_redor2 = plt.subplots()
    for internuc_distance in np.arange(0.4, 0.61, 0.05):
        NTr, S_S0 = redor_curve(internuc_distance, spin_speed, N)
        xvals = NTr * 1000
        ax_redor2.plot(xvals, S_S0, 'k')
        label_curve(ax_redor2, xvals, S_S0, f"{internuc_distance:.2f} nm", idx=25)

    Necos = 20 # for the plot
    ax_redor2.plot(recopl_time[:Necos]*1000, S[:Necos]/S0[:Necos], 'o', alpha=0.8)
    ax_redor2.set_xlabel("Dephasing time (ms)")
    ax_redor2.set_ylabel("S/S0")
    ax_redor2.set_title(f"Señales integradas - ppm range: {ppmRange[0]} to {ppmRange[1]}")
    ax_redor2.set_ylim(0, 1.1)
    ax_redor2.grid(True)
    ax_redor2.legend()

    #=====================================================================
    # Gráficos individuales por par S/S0
    #=====================================================================
    if plot_individual_pairs:
        for i in range(0, len(spec), 2):
            if i + 1 >= len(spec):
                break  # Evita error si hay número impar

            fig_pair, ax_pair = plt.subplots()
            ax_pair.plot(ppmAxis, spec[i, :], label=f"S (index {i})", color='blue')
            ax_pair.plot(ppmAxis, spec[i + 1, :], label=f"S0 (index {i+1})", color='red')
            ax_pair.set_title(f"Par de espectros índice lógico {i//2} - ppm range {ppmRange[0]} to {ppmRange[1]}")
            ax_pair.set_xlabel("ppm")
            ax_pair.set_ylabel("Intensidad")
            ax_pair.set_xlim(np.max(ppmAxis), np.min(ppmAxis))  # orientación NMR
            ax_pair.legend()
            ax_pair.grid(True)
            ax_pair.set_xlim([max(ppmRange), min(ppmRange)])  # Ajuste de límites x
            max_range = np.max(spec[i, np.abs(ppmAxis-0.5*(r1+r2)) < 0.5*(r2-r1)])
            ax_pair.set_ylim([-0.05* max_range, 1.05 * max_range])  # Ajuste de límites y
            ax_pair.axvspan(r1, r2, alpha=0.15, color=color)

# Ajuste final del gráfico general
ax_spec.set_xlim([-75, -90])

# %%
