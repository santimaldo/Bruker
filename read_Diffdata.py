# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate
import matplotlib.ticker as ticker

# info: muestra, expn, ppmRange, bmax
info = ['LiCl-Alcohol-200uL', 12, [-3, 2], 10]
info = ['LiCl-Alcohol-500uL', 1022, [-3, 2], 10]
<<<<<<< Updated upstream
info = ['Li2S6-DME', 3, [-3, 2], 0.6]
info = ['Li2S6-TEGDME', 22, [-3, 2], 100]
# info = ['Li2S6-DME-Delta-10ms', 4, [-1,1], 0.5]
# info = ['Li2S6-DME-Delta-1000ms', 5, [-1,1], 100]
info = ['Li2S6-TEGDME-Delta-100ms', 23, [-3, 2], 100]
=======
info = ['Li2S6-DME-dia1', 3, [-1, 1], 1.5]  # 9/3/23
info = ['Li2S6-DME-dia5', 12, [-0.5, 1], 100]  # 14/3/23

info = ['Li2S6-DME-dia5', 12, [-1, 1], 1]  # 14/3/23


# info = ['Li2S6-TEGDME-dia1', 7, [-2, 2], 100]  # 9/3/23
# info = ['Li2S6-TEGDME-dia5', 15, [-0.5, 0.5], 100]  # 14/3/23
# info = ['Li2S6-TEGDME-dia5', 16, [-0.6, 0.6], 100]  # 14/3/23
>>>>>>> Stashed changes

save = True

# forma del gradiente
gpshape = 'sin'
# factor de correccion: Dref(medido)/Dref(literatura)
factor_b = 1.95
modulo = False
#-------------------- directorios
muestra, expn, ppmRange, bmax = info
# Polisulfuros
path_local = "S:/NMRdata/2022_Polisulfuros/"
path_bruker = f"/2023-03-28_Diff_Polisulfuros/{expn}/"
# Silicio
# path_local = "S:/NMRdata/2022_Silicio/"
# path_bruker = f"2022-12-19_Diff_LiTFSI-SiO2/{expn}/"
# test
# path_local = "S:/NMRdata/2023_tests/"
# path_bruker = f"2023-02-17_Diff_Shimming/{expn}/"

path = path_local + path_bruker
# directorio de guradado
# savepath_local = "G:/Otros ordenadores/Mi PC/"  # Acer
savepath_local = "S:/Posdoc/Li-S/Analisis/2023-03_Li2S6-Diff_DME-TEGDME/"  # Oficina
savepath = f"{savepath_local}"

# --------------------------- Extraigo datos
# 1.7933
datos = DatosProcesadosDiff(path, factor_b=factor_b,
                            bmax=bmax, gpshape=gpshape)
delta = datos.delta
bigDelta = datos.bigDelta
gpmax = datos.gpmax
# -----------------------------------------------


# -----------------------------------------------
msg = f"muestra: {muestra}"
print(msg)
print(f"rango = {info[2]}")
# %%

# --------------------------- grafico un espectro
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

re = datos.espectro.real[1]
im = datos.espectro.imag[1]
if modulo:
    re = np.abs(re+1j*im)

titulo = f"muestra: {muestra}\n" \
         fr"$\Delta = {bigDelta}$ ms, $\delta = {delta}$ ms"

r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
rango = [r1, r2]
fig1d, ax1d = plt.subplots(1, 1)
ax1d.axvspan(r1, r2, alpha=0.2, color='red')
ax1d.axhline(0, color='gray')
ax1d.plot(ppmAxis, re/np.max(re), 'k', lw=2)
ax1d.text(r1-np.abs(0.1*r1), 0.8, "Region de integracion\n(Diff)", color='r')
ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax1d.set_xlabel("NMR Shift [ppm]")
ax1d.set_title(titulo)

plt.axhline(0, color='k')
# -----------------------------------------------

# %%

bvalue, signal = datos.get_Diffdata(ppmRange, absolute=modulo)
bvalue_fit, signal_fit, residuals = datos.Diff1_fit()

smax = np.max(signal)
signal = signal/smax
signal_fit = signal_fit/smax
residuals = residuals/smax

fig, axs = plt.subplots(2, 1, figsize=(6, 7))
# fig.suptitle(muestra)
# -------------
titulo = f"muestra: {muestra}\n" \
         fr"$\Delta = {bigDelta}$ ms, $\delta = {delta}$ ms, "\
         fr"$G_{{max}} = {gpmax:.2f} \%$"
axs[0].set_title(titulo)
axs[0].plot(bvalue, signal, 'ko')
axs[0].plot(bvalue_fit, signal_fit, 'r-')
text = f"$D =$ {datos.fit_results[0]:.4f} $10^{{-9}}m^2/s$ \n " \
       f"$u_D =$ {datos.fit_results[1]:.1g} $10^{{-9}}m^2/s$ \n " \
       f"$r^2 =$ {datos.fit_results[2]:.5g}"
axs[0].text(bvalue[-1]*0.6, 0.5, text,
            multialignment="left")
axs[0].set(xlabel=r'$b_{value} [10^9 s/m^2]$', ylabel=r'$S$')
axs[0].set_yscale('log')
axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1g'))
axs[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1g'))
# -------------
axs[1].plot(bvalue, residuals, 'ko')
axs[1].axhline(0, color='k', linestyle='--')
axs[1].set(xlabel=r'$b_{value} [10^9 s/m^2]$', ylabel=r'Residuos')
axs[1].set_yscale('symlog')
ylim = 10**np.ceil(np.log10(np.max(abs(residuals))))
axs[1].set_ylim([-ylim, ylim])


for ax in axs.flat:
    # ax.label_outer()
    ax.tick_params(direction='in', which='both')


# guardo data:
if save:
    filename0 = f"{muestra}"

    filename = f'{savepath}{filename0}_Diff.png'
    fig.savefig(filename)   # save the figure to file

    filename = f'{savepath}' \
               f'{filename0}_Diff-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n"\
             f"Rango de integracion: {ppmRange} ppm\n"\
             f"bvalue (10^-9 m^2/s)\t S (norm)"
    Diffdata = np.array([bvalue, signal]).T
    np.savetxt(f"{savepath}{filename0}_Diff.dat",
               Diffdata, header=header)

    # data = np.array([ppmAxis, re, im]).T
    # filename = f"{savepath}datos_Diff/primerEspectro/" \
    #            f"{filename0}_primerEspectro.dat"
    # np.savetxt(filename, data)
