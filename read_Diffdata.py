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


# directorio de datos
expn = 6
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_bruker = f"/2022-09-29_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/datos/"
# Nmuestra es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
Nmuestra = 13

expnums = [3, 6]
Nssssss = [16, 13]

save = True
# rango de integracion
ppmRange = [-0.1, 0.2]


# --------------------------- grafico un espectro
datos = DatosProcesadosDiff(path, factor_b=1)
delta = datos.delta
bigDelta = datos.bigDelta

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

re = datos.espectro.real[1]
im = datos.espectro.imag[1]


plt.figure(7532)
plt.plot(ppmAxis, re)
# plt.plot(ppmAxis, im)
plt.xlim(np.max(ppmAxis), np.min(ppmAxis))
plt.axhline(0, color='k')
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
       f"$u_D =$ {datos.fit_results[1]:.1g} \n " \
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

    header = f"Archivo de datos: {path_bruker}\n"\
             f"Rango de integracion: {ppmRange} ppm\n"\
             f"bvalue (10^-9 m^2/s)\t S (norm)"
    Diffdata = np.array([bvalue, signal]).T
    np.savetxt(f"{savepath}/{filename0}_Diff.dat", Diffdata, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{filename0}_primerEspectro.dat", data)
