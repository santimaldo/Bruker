# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

from paths import NMRDATA, HOME
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


def fit_nutation(x, signal, n_slices):
    """
    Fit a nutation curve with a sine.

    Parameters
    ----------
    x : ndarray
        Pulse lengths (p1) in us.
    signal : ndarray
        Real-valued nutation curve (integrals vs p1).
    n_slices : int
        Number of slices, used for the initial frequency guess.

    Returns
    -------
    popt : ndarray
        Fitted [A, f, phi, C] (f in MHz), or NaNs if the fit fails.
    """
    A0 = (np.max(signal)-np.min(signal))/2
    C0 = np.mean(signal)
    f0 = 1/(2*n_slices)
    phi0 = 0
    try:
        popt, pcov = curve_fit(sine, x, signal, p0=[A0, f0, phi0, C0],
                               maxfev=100000)
    except RuntimeError:
        popt = np.full(4, np.nan)
    return popt

################## end Functions #######################


absolute = False
autophase = False

path = rf"{NMRDATA}\400dnp\2026-09-10_3.2mm_LiMetal_LiF_RC1/"
savepath = rf"{HOME}\LiMetal\DNP\analysis\2026-09_Limetal_LiF\10kHz/"

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
p1 = np.arange(0.25, 80.1, 0.25)
# expn = 226999
# p1 = np.arange(0.25, 5.1, 0.25)
# BSMS del experimento particular.
# Se usa UNICAMENTE para calibrar el eje ppm.
bsms = -3500
muestra = "0kHz_nutation"
save = False
print("WARNING!!!!!!!!!! EL P1 ESTA HARDOCODEADO")
plotRange = [500, -200]
# Rango general usado para definir los bins de ppm y Ancho de cada intervalo de integracion
# ppmRange, ppm_bin_width = [301, 0], 50
#ppm_bin_width, ppmRanges = None, [[275, 240], [180, 50]]
ppm_bin_width, ppmRanges = None, [[1200,-800]]#[[300,200],[200,0]]



######### 6Li
path = rf"{NMRDATA}\400dnp\2026-10-22_3.2mm_LiMetal_LiF_RC3/"
savepath = rf""

expn = 308998
p1 = np.arange(1, 64.1, 1)
bsms = 4250
muestra = "0kHz_nutation"
save = False
print("WARNING!!!!!!!!!! EL P1 ESTA HARDOCODEADO")
plotRange = [400, -200]
# Rango general usado para definir los bins de ppm y Ancho de cada intervalo de integracion
# ppmRange, ppm_bin_width = [301, 0], 50
#ppm_bin_width, ppmRanges = None, [[275, 240], [180, 50]]
ppm_bin_width, ppmRanges = None, [[300,200],[200,50]]



# expn = 300999
# p1 = np.arange(2, 32.1, 2)
# bsms = 0
# ppm_bin_width, ppmRanges = None, [[280,240]]




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

# Parte real: integral del espectro real; parte imaginaria: integral del
# espectro imaginario
integrals = np.zeros((n_slices, len(ppmRanges))).astype("complex")


for kk in range(n_slices):

    spec1d = spec[kk, :]
    speci1d = speci[kk, :]

    if absolute:
        spec1d = np.abs(spec1d + 1j*speci1d)
        speci1d = np.zeros_like(speci1d)  # Im: cero en modo absoluto

    elif autophase:
        spec_phased = ng.proc_autophase.autops(spec1d+1j*speci1d, "acme")
        spec1d = spec_phased.real
        speci1d = spec_phased.imag  # Im: espectro imaginario fasado

    # Integracion en cada intervalo de ppm
    for ii, ppmRange_i in enumerate(ppmRanges):

        r1, r2 = [np.min(ppmRange_i), np.max(ppmRange_i)]

        condition = np.abs(ppmAxis-np.mean(ppmRange_i)) < np.abs(r1-r2)/2

        spec_to_integrate = (spec1d+1j*speci1d)[condition]  # Re + i*Im
        ppm_to_integrate = ppmAxis[condition]

        # Los doy vuelta porque el eje ppm esta descendiendo
        ppm_to_integrate = ppm_to_integrate[::-1]
        spec_to_integrate = spec_to_integrate[::-1]

        # Integral
        integrals[kk, ii] = simpson(spec_to_integrate, x=ppm_to_integrate)  # Re: integral real; Im: integral imaginaria


#=====================================================================
# Ajuste de las curvas de nutacion
#=====================================================================

frequencies = np.zeros(len(ppmRanges))
fit_params = np.zeros((len(ppmRanges), 4))
frequencies_imag = np.zeros(len(ppmRanges))  # Im: frecuencias del ajuste
fit_params_imag = np.zeros((len(ppmRanges), 4))  # Im: parametros del ajuste

for ii, ppmRange_i in enumerate(ppmRanges):

    fit_params[ii, :] = fit_nutation(p1, np.real(integrals[:, ii]), n_slices)
    frequencies[ii] = abs(fit_params[ii, 1])

    fit_params_imag[ii, :] = fit_nutation(p1, np.imag(integrals[:, ii]), n_slices)  # Im: ajuste
    frequencies_imag[ii] = abs(fit_params_imag[ii, 1])  # Im: frecuencia


#=====================================================================
# Graficar las curvas de nutacion y sus ajustes
#=====================================================================

cmap = plt.get_cmap("gnuplot")
norm_ppmRange = mpl.colors.Normalize(vmin=0, vmax=len(ppmRanges)-1)

fig, ax = plt.subplots(figsize=(6, 8))

for ii, ppmRange_i in enumerate(ppmRanges):

    color = cmap(norm_ppmRange(ii))
    signal = np.real(integrals[:, ii])
    signal_imag = np.imag(integrals[:, ii])  # Im: curva de nutacion
    ppm_center = np.mean(ppmRange_i)

    ax.plot(p1, signal, 'o', color=color, label=f"Re: {ppm_center:.0f} ppm")
    ax.plot(p1, signal_imag, 'o', color=color, markerfacecolor="white",  # Im: datos
            label=f"Im: {ppm_center:.0f} ppm")

    if not np.isnan(frequencies[ii]):
        ax.plot(p1, sine(p1, *fit_params[ii, :]), '-', color=color)
    if not np.isnan(frequencies_imag[ii]):  # Im: ajuste
        ax.plot(p1, sine(p1, *fit_params_imag[ii, :]), ':', color=color)  # Im: ajuste

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

    ax_i.plot(p1, np.real(integrals[:, ii]), 'o', color="tab:blue", label="Re")
    ax_i.plot(p1, np.imag(integrals[:, ii]), 'o', color="tab:red",  # Im: datos
              markerfacecolor="white", label="Im")

    if not np.isnan(frequencies[ii]):
        ax_i.plot(p1, sine(p1, *fit_params[ii, :]), '-', color="tab:blue")
    if not np.isnan(frequencies_imag[ii]):  # Im: ajuste
        ax_i.plot(p1, sine(p1, *fit_params_imag[ii, :]), '-', color="tab:red")  # Im: ajuste
    ax_i.legend(fontsize=8)

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

ax_freq.plot(ppm_centers, frequencies*1e3, 'o-', label="Re")
ax_freq.plot(ppm_centers, frequencies_imag*1e3, 'o--', markerfacecolor="white",  # Im: frecuencia
             label="Im")
ax_freq.legend()

ax_freq.set_xlabel(r"$\delta$ (ppm)")
ax_freq.set_ylabel("Nutation frequency [kHz]")
ax_freq.set_title(f"Nutation frequency - expn {expn}")
ax_freq.set_xlim(max(ppmRange_i), min(ppmRange_i))
ax_freq.grid(alpha=0.2)


#=====================================================================
# Espectro de referencia
#=====================================================================

fig_spec, ax_spec = plt.subplots(figsize=(8, 5))
ax_spec.axhline(0, color="k")
ax_spec.plot(ppmAxis, spec[4, :], label="slice 4 (Re)")
ax_spec.plot(ppmAxis, speci[4, :], label="slice 4 (Im)")  # Im: espectro
ax_spec.legend()

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


#%%=====================================================================
# Transformada de Fourier de las curvas de nutacion
#=====================================================================

# p1 en us -> rfftfreq da MHz -> *1e3 para kHz
zero_fill = 4  # factor de zero-filling (1 = sin zero-filling)
n_fft = zero_fill*len(p1)
nut_freqs_kHz = np.fft.rfftfreq(n_fft, d=p1[1]-p1[0]) * 1e3

nut_fft = np.zeros((len(ppmRanges), len(nut_freqs_kHz)))
nut_fft_imag = np.zeros((len(ppmRanges), len(nut_freqs_kHz)))  # Im: FFT
for ii, ppmRange_i in enumerate(ppmRanges):
    signal = np.real(integrals[:, ii])
    signal_imag = np.imag(integrals[:, ii])  # Im: curva de nutacion
    # Restar el valor medio para eliminar la componente DC
    signal = signal - np.mean(signal)
    signal_imag = signal_imag - np.mean(signal_imag)  # Im: sin componente DC
    nut_fft[ii, :] = np.abs(np.fft.rfft(signal, n=n_fft))
    nut_fft_imag[ii, :] = np.abs(np.fft.rfft(signal_imag, n=n_fft))  # Im: FFT

fig_fft, ax_fft = plt.subplots(figsize=(6, 4))
for ii, ppmRange_i in enumerate(ppmRanges):
    color = cmap(norm_ppmRange(ii))
    # ambas normalizadas al maximo de la parte real, para comparar amplitudes
    norm_fft = np.max(nut_fft[ii, :])
    ax_fft.plot(nut_freqs_kHz, nut_fft[ii, :]/norm_fft, '-',
                color=color,
                label=f"Re: {ppmRange_i[0]:.0f} - {ppmRange_i[1]:.0f} ppm")
    ax_fft.plot(nut_freqs_kHz, nut_fft_imag[ii, :]/norm_fft, ':',  # Im: FFT
                color=color,
                label=f"Im: {ppmRange_i[0]:.0f} - {ppmRange_i[1]:.0f} ppm")
    # frecuencias del ajuste, para comparar
    if not np.isnan(frequencies[ii]):
        ax_fft.axvline(frequencies[ii]*1e3, color=color, ls="--", lw=1)
    if not np.isnan(frequencies_imag[ii]):  # Im: frecuencia del ajuste
        ax_fft.axvline(frequencies_imag[ii]*1e3, color=color, ls="-.", lw=1)  # Im: frecuencia del ajuste

ax_fft.set_xlabel("Nutation frequency [kHz]")
ax_fft.set_ylabel("|FFT| (normalized)")
ax_fft.set_title(f"FFT of the nutation curves - expn {expn}")
ax_fft.legend(fontsize=8)
ax_fft.grid(alpha=0.2)
ax_fft.set_xlim(0,400)
fig_fft.tight_layout()


# %%
