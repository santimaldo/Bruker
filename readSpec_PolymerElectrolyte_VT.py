# -*- coding: utf-8 -*-
"""
Created on Jul 29 2025

@author: Santi

Extrae y grafica espectros Bruker adquiridos a diferentes temperaturas
(Variable Temperature experiment) usando un colormap tipo plasma.
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

# === Configuración general ===
nucleo = "1H"
nexp_base = 36
T_nominal = np.arange(-20, 41, 10)
T_real = np.array([-17.8, -9.9, 0, 10, 20, 30, 40])
usar_Treal = True  # <--- Cambiar a False para usar T_nominal

path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-07-16_PEO-solid-electrolyte_VT/"
save = False
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2025-07_VT/"
savepath_especifico = ""

ppmRange = [20, -20]

# === Generación de experimentos y temperaturas ===
expns = [nexp_base + 10 * i for i in range(len(T_nominal))]
temperaturas = T_real if usar_Treal else T_nominal

# === Preparar colores y figura ===
cmap = plt.get_cmap("jet_r")
colores = [cmap(i * 0.9 / (len(temperaturas) - 1)) for i in range(len(temperaturas))]

fig, ax = plt.subplots(num=1, nrows=1, ncols=1)

for i, (expn, temp) in enumerate(zip(np.flip(expns), np.flip(temperaturas))):
    datos = DatosProcesados(f'{path}{expn}/')

    if ppmRange is not None:
        datos.espectro.ppmSelect(ppmRange)

    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis

    color = colores[i]
    ax.plot(ppmAxis, re, linewidth=2, color=color, label=f"{temp} °C")

    if save:
        savepath = f"{savepath_local}{savepath_especifico}"
        header = "ppmAxis\t real \t imag"
        dataexport = np.array([ppmAxis, re, im]).T
        filename = f'{savepath}/{nucleo}_VTexp{expn}_{temp}C.dat'
        np.savetxt(filename, dataexport, header=header)

# === Formato del gráfico ===
ax.axhline(0, color='k', ls='--')
ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
ax.set_ylabel("Intensidad (u.a.)")
ax.set_title("Variable Temperature NMR")
ax.invert_xaxis()
ax.legend(title="Temperature")

plt.tight_layout()
plt.show()
