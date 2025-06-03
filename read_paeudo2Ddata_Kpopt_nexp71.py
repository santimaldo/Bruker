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
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

################## Functions ###########################

def read_bsms_field(path_archivo):
    with open(path_archivo+'Klog_opt', 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                x = eval(contenido)
                return np.array(x)
    raise ValueError("No se encontro una linea que empiece con 'x='.")

def integrate_two_windows_from_minimum(ppmAxis, spec1d, width_ppm):
    smoothed = gaussian_filter1d(spec1d, sigma=3)
    peaks, _ = find_peaks(smoothed)
    if len(peaks) < 2:
        raise ValueError("No se encontraron dos máximos para definir las ventanas.")
    peak_vals = smoothed[peaks]
    idx_peaks = peaks[np.argsort(peak_vals)[-2:]]
    idx_peaks.sort()
    idx_min = np.argmin(smoothed[idx_peaks[0]:idx_peaks[1]+1]) + idx_peaks[0]
    ppm_min = ppmAxis[idx_min]
    r1 = (ppm_min - width_ppm, ppm_min)
    r2 = (ppm_min, ppm_min + width_ppm)
    mask1 = (ppmAxis >= min(r1)) & (ppmAxis <= max(r1))
    mask2 = (ppmAxis >= min(r2)) & (ppmAxis <= max(r2))
    I1 = -integrate.simpson(spec1d[mask1], x=ppmAxis[mask1])
    I2 = -integrate.simpson(spec1d[mask2], x=ppmAxis[mask2])
    return I1, I2, r1, r2

################## end Functions #######################

# parámetros
expn = 71
absolute = False
path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\InSitu_May_2025/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker/analysis/2025-05_mesh/spec1d_MeshAndLFP/"
muestra = ""
save = False
ppmRange = [200, -200]
plotRange = [100, -100]

window_width = 50  # ancho de cada ventana desde el mínimo local

# iniciar figuras
fig_spec, ax_spec = plt.subplots(num=17856)
fig_popt, ax_popt = plt.subplots(num=382910)

# leer datos
path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(path_2D, read_pp=False)
datos.espectro.ppmSelect(ppmRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real
speci = datos.espectro.imag
bsms_field = read_bsms_field(path_2D)

integrals_1 = []
integrals_2 = []

for kk in range(bsms_field.size):
    spec1d = spec[kk, :]
    speci1d = speci[kk, :]

    if absolute:
        abs_spec = np.abs(spec1d + 1j * speci1d)
        spec1d = abs_spec
        spec[kk, :] = abs_spec

    ax_spec.plot(ppmAxis, spec1d)

    I1, I2, r1, r2 = integrate_two_windows_from_minimum(ppmAxis, spec1d, width_ppm=window_width)
    integrals_1.append(I1)
    integrals_2.append(I2)
    ax_spec.axvspan(*r1, alpha=0.15, color='tab:blue')
    ax_spec.axvspan(*r2, alpha=0.15, color='tab:orange')
ax_spec.set_xlim(plotRange)

Signals = np.array([integrals_1, integrals_2]).T

if save:
    for kk in range(Signals.shape[0]):
        np.savetxt(f"{savepath}/{muestra}_bsms_{bsms_field[kk]}.dat",
                np.array([ppmAxis, spec[kk, :]]).T,
                header="ppmAxis\treal")

ax_popt.plot(bsms_field, integrals_1, 'o', label="Peak 1")
ax_popt.plot(bsms_field, integrals_2, 's', label="Peak 2")
ax_popt.legend()

np.savetxt(f"{savepath}/{muestra}_integrals_vs_bsms.dat",
           np.column_stack([bsms_field, integrals_1, integrals_2]),
           header="bsms_field\tIntegral1\tIntegral2")

print(np.array([bsms_field, integrals_1, integrals_2]).T)
