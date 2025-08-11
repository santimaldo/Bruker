# -*- coding: utf-8 -*-
"""
Created on Jul 29 2025

@author: Santi

Extrae y grafica espectros Bruker adquiridos a diferentes temperaturas
(Variable Temperature experiment)
"""

import nmrglue as ng
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import numpy as np
from Datos import *

# === Configuración general ===
nucleo = "19F"
nexp_base = 27
T_nominal = np.array([-19, -10, 0, 10, 20, 30, 40])

path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-08-05_PEO-PTT-solid-electrolyte_VT/"
save = True
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2025-08_500MHz_VT/"
savepath_especifico = "LiTFSI-PEO-PTT"

ppmRange = [-10, -150]

# === Generación de experimentos y temperaturas ===
expns = [nexp_base + 10 * i for i in range(len(T_nominal))]
temperaturas = T_nominal

# === Preparar colores y figura ===
cmap = plt.get_cmap("jet_r")
colores = [cmap(i * 0.9 / (len(temperaturas) - 1)) for i in range(len(temperaturas))]

fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=(4,3.25))

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
ax.set_xlabel(f"{nucleo} Chemical Shift [ppm]")
ax.set_ylabel("Intensity (a.u.)")
# ax.set_title("Variable Temperature NMR")
ax.invert_xaxis()
ax.set_yticks([])
# ax.legend(title="Temperature")
maximo = 1e3 
ax.set_ylim(-0.05*maximo, 1.05*maximo)

plt.tight_layout()
plt.show()
