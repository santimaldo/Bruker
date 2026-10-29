# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

from paths import NMRDATA_LOCAL, PROJECTS
import nmrglue as ng
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from Datos import *
from scipy.integrate import simpson
from scipy.optimize import curve_fit

################## Functions ###########################

def sine(x, A, f, phi, C):
    return A*np.sin(2*np.pi*f*x + phi) + C

################## end Functions #######################


absolute = False
autophase = False

path = rf"{NMRDATA_LOCAL}\400dnp\2026-09-10_3.2mm_LiMetal_LiF_RC1/"
savepath = rf"{PROJECTS}\LiMetal\DNP\analysis\2026-09_Limetal_LiF\10kHz/"

# expn = 41999
# # BSMS del experimento particular.
# # Se usa UNICAMENTE para calibrar el eje ppm.
# bsms = -3750
# muestra = "10kHz_nutation"
# save = False
# p1 = np.arange(0.1, 7.91, 0.15)
# print("WARNING!!!!!!!!!! EL P1 ESTA HARDOCODEADO")
# plotRange = [300, -20]
# # Rango general usado para definir los bins de ppm
# ppmRange = [260, 0]
# # Ancho de cada intervalo de integracion
# ppm_bin_width = 50

expn = 260999
# BSMS del experimento particular.
# Se usa UNICAMENTE para calibrar el eje ppm.
bsms = -3500
muestra = "0kHz_nutation"
save = False
p1 = np.arange(0.25, 80.1, 0.25)
print("WARNING!!!!!!!!!! EL P1 ESTA HARDOCODEADO")
plotRange = [300, -20]
# Rango general usado para definir los bins de ppm y Ancho de cada intervalo de integracion
# ppmRange, ppm_bin_width = [301, 991], 100
ppm_bin_width, ppmRanges = None, [[275, 240], [180, 50]]

#=====================================================================
# Leer experimento 2D
#=====================================================================

path_2D = f"{path}/{expn}/"

datos = DatosProcesados2D(f'{path}/{expn}/', read_pp=False)

spec = datos.espectro.real
speci = datos.espectro.imag

# Cada fila del experimento 2D corresponde a un indice de p1.
# NO se lee nada relacionado con BSMS.
n_slices = datos.acqu2s.TD


#=====================================================================
# Calibrar el eje ppm
#=====================================================================

# La calibracion del eje ppm se mantiene igual que en el codigo original,
# pero ahora usando el BSMS ingresado manualmente para este experimento.
ppmAxis = datos.espectro.ppmAxis - (-1.1665e-2) * bsms


#=====================================================================
# Definir los intervalos de ppm
#=====================================================================
if ppm_bin_width is not None:
    ppm_bins = np.arange(ppmRange[0], ppmRange[1], -ppm_bin_width)
    ppmRanges = [[ppm+ppm_bin_width/2, ppm-ppm_bin_width/2] for ppm in ppm_bins]


#=====================================================================
# Integracion para cada ppm y cada indice de p1
#=====================================================================

integrals = np.zeros((n_slices, len(ppmRanges))).astype("complex")

for kk in range(n_slices):

    spec1d = spec[kk, :]
    speci1d = speci[kk, :]

    if absolute:
        spec1d = np.abs(spec1d + 1j*speci1d)

    elif autophase:
        spec1d = ng.proc_autophase.autops(spec1d+1j*speci1d, "acme").real

    # Integracion en cada intervalo de ppm
    for ii, ppmRange_i in enumerate(ppmRanges):

        r1, r2 = [np.min(ppmRange_i), np.max(ppmRange_i)]

        condition = np.abs(ppmAxis-np.mean(ppmRange_i)) < np.abs(r1-r2)/2

        spec_to_integrate = (spec1d+1j*speci1d)[condition]
        ppm_to_integrate = ppmAxis[condition]

        # Los doy vuelta porque el eje ppm esta descendiendo
        ppm_to_integrate = ppm_to_integrate[::-1]
        spec_to_integrate = spec_to_integrate[::-1]

        # Integral
        integrals[kk, ii] = simpson(spec_to_integrate, x=ppm_to_integrate)


#=====================================================================
# Ajuste de las curvas de nutacion
#=====================================================================

frequencies = np.zeros(len(ppmRanges))
fit_params = np.zeros((len(ppmRanges), 4))

for ii, ppmRange_i in enumerate(ppmRanges):

    signal = np.real(integrals[:, ii])

    # Estimacion inicial
    A0 = (np.max(signal)-np.min(signal))/2
    C0 = np.mean(signal)
    f0 = 1/(2*n_slices)
    phi0 = 0

    p0 = [A0, f0, phi0, C0]

    try:
        popt, pcov = curve_fit(sine, p1, signal, p0=p0, maxfev=100000)
        fit_params[ii, :] = popt
        frequencies[ii] = abs(popt[1])

    except RuntimeError:
        fit_params[ii, :] = np.nan
        frequencies[ii] = np.nan


#=====================================================================
# Graficar las curvas de nutacion y sus ajustes
#=====================================================================

cmap = plt.get_cmap("gnuplot")
norm_ppmRange = mpl.colors.Normalize(vmin=0, vmax=len(ppmRanges)-1)

fig, ax = plt.subplots(figsize=(6, 8))

for ii, ppmRange_i in enumerate(ppmRanges):

    color = cmap(norm_ppmRange(ii))
    signal = np.real(integrals[:, ii])
    ppm_center = np.mean(ppmRange_i)

    ax.plot(p1, signal, 'o', color=color, label=f"{ppm_center:.0f} ppm")

    if not np.isnan(frequencies[ii]):
        ax.plot(p1, sine(p1, *fit_params[ii, :]), '-', color=color)

ax.set_xlabel("p1 [us]")
ax.set_ylabel("Integral [a.u.]")
ax.set_title(f"Curvas de nutación - expn {expn}")
ax.legend(title="ppm")
ax.grid(alpha=0.2)


#=====================================================================
# Grafico individual para cada intervalo de ppm
#=====================================================================

fig_ranges, axs_ranges = plt.subplots(nrows=len(ppmRanges), figsize=(6, 2*len(ppmRanges)), squeeze=False)

for ii, ppmRange_i in enumerate(ppmRanges):

    ax_i = axs_ranges[ii, 0]
    ppm_center = np.mean(ppmRange_i)

    ax_i.plot(p1, np.real(integrals[:, ii]), 'o')

    if not np.isnan(frequencies[ii]):
        ax_i.plot(p1, sine(p1, *fit_params[ii, :]), '-')

    ax_i.set_ylabel("Integral")
    ax_i.set_title(f"{ppmRange_i[0]:.0f} - {ppmRange_i[1]:.0f} ppm")
    ax_i.grid(alpha=0.2)

axs_ranges[-1, 0].set_xlabel("p1 index")
fig_ranges.tight_layout()


#=====================================================================
# Frecuencia de nutacion vs ppm
#=====================================================================

ppm_centers = np.array([np.mean(ppmRange_i) for ppmRange_i in ppmRanges])

fig_freq, ax_freq = plt.subplots()

ax_freq.plot(ppm_centers, frequencies*1e3, 'o-')

ax_freq.set_xlabel(r"$\delta$ (ppm)")
ax_freq.set_ylabel("Nutation frequency [kHz]")
ax_freq.set_title(f"Nutation frequency - expn {expn}")
ax_freq.set_xlim(max(ppmRange), min(ppmRange))
ax_freq.grid(alpha=0.2)


#=====================================================================
# Espectro de referencia
#=====================================================================

fig_spec, ax_spec = plt.subplots(figsize=(8, 5))

ax_spec.plot(ppmAxis, spec[0, :], label="slice 0")

ax_spec.set_xlim(max(plotRange), min(plotRange))
ax_spec.set_xlabel(r"$\delta$ (ppm)")
ax_spec.set_ylabel("Signal [a.u.]")
ax_spec.set_title(f"Espectro - expn {expn}")

for ppmRange_i in ppmRanges:

    r1, r2 = [np.min(ppmRange_i), np.max(ppmRange_i)]
    ax_spec.axvspan(r1, r2, alpha=0.1)




######################################################################
######################################################################
######################################################################
######################################################################
# #%%=====================================================================
# # Transformada de Fourier de las curvas de nutacion
# #=====================================================================

# fft_spectra = np.zeros((len(ppmRanges), len(p1)//2+1))
# fft_freqs = np.fft.rfftfreq(len(p1), d=p1[1]-p1[0])


# for ii, ppmRange_i in enumerate(ppmRanges):

#     signal = np.real(integrals[:, ii])

#     # Restar el valor medio para eliminar la componente DC
#     signal = signal - np.mean(signal)

#     # Transformada de Fourier
#     fft_spectra[ii, :] = np.abs(np.fft.rfft(signal))

# #=====================================================================
# # Grafico waterfall de las transformadas de Fourier
# #=====================================================================

# fig_fft, ax_fft = plt.subplots(figsize=(8, 8))

# offset = 0

# for ii, ppmRange_i in enumerate(ppmRanges):

#     ppm_center = np.mean(ppmRange_i)

#     ax_fft.plot(fft_freqs, fft_spectra[ii, :]/max(fft_spectra[ii, :]),
#                 label=f"{ppm_center:.0f} ppm")

# # ax_fft.set_xscale("log")
# ax_fft.set_xlabel("Frequency [MHz]")
# ax_fft.set_ylabel("FFT intensity + offset [a.u.]")
# ax_fft.set_title(f"FFT de las curvas de nutacion - expn {expn}")
# ax_fft.grid(alpha=0.2)


# %%
