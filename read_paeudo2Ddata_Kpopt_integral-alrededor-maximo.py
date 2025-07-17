# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Autophase import autophase
import scipy.integrate as integrate
import re

################## Functions ###########################

def read_bsms_field(path_archivo):
    """
    Lee un archivo con formato específico y devuelve el array 'x' como lista de floats.
    """
    with open(path_archivo + 'Klog_opt', 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                x = eval(contenido)
                return x
    raise ValueError("No se encontró una línea que empiece con 'x='.")

################## end Functions #######################

# Configuración
expns = np.arange(130, 135)
absolute = False
autoph = True
path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/2025-06-27_InSitu/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\Li-on-mesh-LFP"

nucleo = "7Li"
muestra = ""
save = False
plotRange = [600, -100]
window_width = 50  # ancho de integración centrado en el máximo (en ppm)

# Gráficos
fig_spec, ax_spec = plt.subplots(num=17856)
fig_popt, ax_popt = plt.subplots(num=382910)

for expn in expns:
    Signals = np.array([])

    # Lectura y procesamiento de datos
    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados2D(path_2D, read_pp=False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    speci = datos.espectro.imag

    for kk in range(spec.shape[0]):
        spec1d = spec[kk, :]
        speci1d = speci[kk, :]

        if absolute:
            abs_spec = np.abs(spec1d + 1j * speci1d)
            spec[kk, :] = abs_spec
            spec1d = abs_spec

        elif autoph:
            spec1d = ng.proc_autophase.autops(spec1d + 1j * speci1d, "acme")
            spec[kk, :] = ng.process.proc_bl.cbf(spec1d.real, last=100)

        ax_spec.plot(ppmAxis, spec[kk, :])

    bsms_field = read_bsms_field(path_2D)

    # Integración centrada en el máximo
    signal = np.zeros(spec.shape[0])
    for kk in range(spec.shape[0]):
        spectrum_kk = spec[kk, :]
        max_index = np.argmax(np.abs(spectrum_kk))
        center_ppm = ppmAxis[max_index]
        r1 = center_ppm - window_width / 2
        r2 = center_ppm + window_width / 2

        # Adaptar a eje decreciente
        lower, upper = sorted([r1, r2], reverse=True)
        mask = (ppmAxis <= lower) & (ppmAxis >= upper)

        signal[kk] = -integrate.simpson(spectrum_kk[mask], x=ppmAxis[mask])

        # Marcar solo una vez
        if kk == 0 and expn == expns[0]:
            ax_spec.axvspan(lower, upper, alpha=0.3, color='orange')
            ax_spec.axvline(center_ppm, linestyle='--', color='darkred', alpha=0.7)

    Signals = signal[signal != 0]  # Eliminar ceros si se detuvo la adquisición

    if expn == expns[0]:
        all_Signals = Signals
        all_bsms_field = bsms_field
    else:
        all_Signals = np.concatenate((all_Signals, Signals))
        all_bsms_field = np.concatenate((all_bsms_field, bsms_field))

    if save:
        for kk in range(Signals.size):
            np.savetxt(f"{savepath}/{muestra}_bsms_{bsms_field[kk]}.dat",
                       np.array([ppmAxis, spec[kk, :]]).T,
                       header="ppmAxis\treal")

    ax_popt.plot(bsms_field, Signals, 'o')

# Resultado final
print(np.array([all_bsms_field, all_Signals]).T)
