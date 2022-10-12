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

# Nmuestra es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
Nmuestra = 18
# directorio de datos:
fecha = '10-03'  # 'MM-DD'
expn = 14
# rango de integracion
ppmRange = [-1, 1.5]
save = False
bmax = np.inf


#-------------------- directorios
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_bruker = f"/2022-{fecha}_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/datos/"


# --------------------------- Extraigo datos
datos = DatosProcesadosDiff(path, factor_b=1, bmax=bmax)
delta = datos.delta
bigDelta = datos.bigDelta
# -----------------------------------------------


# -----------------------------------------------
# datos de la muestra
muestra = f"M{Nmuestra}"
N = Nmuestra-10
# matriz es el medio: Bulk, Q30 o Q3
if N < 3:
    matriz = 'Q3'
elif N < 6:
    matriz = 'Q30'
elif N < 9:
    matriz = 'Bulk'
# pc es el porcentaje de glicerol: 50%, 70% o 90%
if N % 3 == 0:
    pc = 50
elif N % 3 == 1:
    pc = 70
elif N % 3 == 2:
    pc = 90
msg = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}"
print(msg)
# %%

# --------------------------- grafico un espectro
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

re = datos.espectro.real[1]
im = datos.espectro.imag[1]

titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n" \
         fr"$\Delta = {bigDelta}$ ms, $\delta = {delta}$ ms"

r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
fig1d, ax1d = plt.subplots(1, 1, num=7532)
ax1d.axvspan(r1, r2, alpha=0.2, color='red')
ax1d.axhline(0, color='gray')
ax1d.plot(ppmAxis, re/np.max(re), 'k', lw=2)
ax1d.text(r1-1, 0.8, "Region de integracion (Diff)", color='r')
ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax1d.set_xlabel("NMR Shift [ppm]")
ax1d.set_title(titulo)

plt.axhline(0, color='k')
# -----------------------------------------------

# %%

bvalue, signal = datos.get_Diffdata(ppmRange)
bvalue_fit, signal_fit, residuals = datos.Diff1_fit()

smax = np.max(signal)
signal = signal/smax
signal_fit = signal_fit/smax
residuals = residuals/smax

fig, axs = plt.subplots(2, 1, figsize=(6, 7))
# fig.suptitle(muestra)
# -------------
titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n" \
         fr"$\Delta = {bigDelta}$ ms, $\delta = {delta}$ ms"
axs[0].set_title(titulo)
axs[0].plot(bvalue, signal, 'ko')
axs[0].plot(bvalue_fit, signal_fit, 'r-')
text = f"$D =$ {datos.fit_results[0]:.3f} $10^{{-9}}m^2/s$ \n " \
       f"$u_D =$ {datos.fit_results[1]:.1g} $10^{{-9}}m^2/s$ \n " \
       f"$r^2 =$ {datos.fit_results[2]:.5g}"
axs[0].text(bvalue[-1]*0.6, 0.7*signal[0], text,
            multialignment="left")
axs[0].set(xlabel=r'$b_{value} [10^9 s^2/m]$', ylabel=r'$S$')
axs[0].set_yscale('log')
axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1g'))
axs[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1g'))
# -------------
axs[1].plot(bvalue, residuals, 'ko')
axs[1].axhline(0, color='k', linestyle='--')
axs[1].set(xlabel=r'$b_{value} [10^9 s^2/m]$', ylabel=r'Residuos')
axs[1].set_yscale('symlog')
ylim = 10**np.ceil(np.log10(np.max(abs(residuals))))
axs[1].set_ylim([-ylim, ylim])


for ax in axs.flat:
    # ax.label_outer()
    ax.tick_params(direction='in', which='both')


# guardo data:
if save:
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}/{filename0}_Diff.png'
    fig.savefig(filename)   # save the figure to file

    filename = f'{savepath}/{filename0}_Diff-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n"\
             f"Rango de integracion: {ppmRange} ppm\n"\
             f"bvalue (10^-9 m^2/s)\t S (norm)"
    Diffdata = np.array([bvalue, signal]).T
    np.savetxt(f"{savepath}/{filename0}_Diff.dat", Diffdata, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{filename0}_primerEspectro.dat", data)
