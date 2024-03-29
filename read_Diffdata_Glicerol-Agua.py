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


# Nmuestra:  es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
# fecha: correspondiente al directorio de datos.  'MM-DD'
# expn: numero de archivo
# ppmRange: rango de integracion
# bmax: maximo valor de bvalue utilizado para ajustar
#
# info = [Nmuestra, fecha, expn, ppmRange, bmax]
# Bulk:
info = [999, '2023-03-08', 3, [-0.5,0.5], 0.7]  
# info = [24, '2022-11-14', 4, [-0.5, 0.5], 1.6]
# info = [16, '09-29', 3, [-0.5, 0.5], 3]
# info = [17, '09-29', 10, [-0.5, 0.5], np.inf]
# info = [18, '10-03', 14, [-1, 1], 50]
# Q30:
# info = [22, '11-14', 12, [-2.5, 2.5], 2.5]
# info = [13, '09-29', 6, [-1.3, 0.7], 20]
# info = [13, '11-17', 122, [-1.3,0.7], 100] # re-medicion
# info = [14, '10-03', 3, [-2.5, 2.5], 14]
# info = [15, '10-11', 3, [-2.5, 2.5], np.inf]
# Q3
# info = [21, '11-14', 23, [-2, 2], 200]
# info = [11, '10-20', 20, [-1, 1], np.inf]

# referencia LiCl 0.25 M
# Gmax = 90%:
# info = [0, '11-17', 4, [-2,2], 0.2] # bigDelta= 5ms, delta=0.5ms # nooooooo
# info = [0, '11-17', 9, [-2,2], 0.3] # bigDelta= 5ms, delta=0.5ms # nooooooo
# info = [0, '11-17', 10, [-2,2], 1.25] # bigDelta=20ms, delta=0.2ms # no
# info = [0, '11-17', 12, [-2,2], 10] # bigDelta=20ms, delta=0.36ms # no
# info = [0, '11-17', 14, [-1,1], 10] # bigDelta=20ms, delta=0.36ms
# Gmax < 16%:
# info = [0, '11-17', 3, [-0.5, 0.5], 200]
# factor_b = 1.88  # bigDelta=10ms, delta=2.0ms
# info = [0, '11-17', 3, [-0.5, 0.5], 200]  # bigDelta=10ms, delta=2.0ms
# info = [0, '11-17', 6, [-0.5, 0.5], 200]  # bigDelta=10ms, delta=4.0ms
# info = [0, '11-17', 5, [-0.5, 0.5], 200]  # bigDelta=50ms, delta=2.0ms
# info = [0, '11-17', 7, [-0.5, 0.5], 200]  # bigDelta=50ms, delta=4.0ms
# Juguito de la Q13
# info = [0, '11-17', 32, [-1,2],200]

save = True
Nmuestra, fecha, expn, ppmRange, bmax = info
#-------------------- directorios
# path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_local = "S:/NMRdata/2022_Glicerol-Agua_CNEA/"
path_bruker = f"{fecha}_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath_local = "G:/Otros ordenadores/Oficina/" # Acer
# savepath_local = "S:/"  # Oficina
savepath = f"{savepath_local}Posdoc/CNEA/Glicerol-Agua/analisis/7Li/"

gpshape = 'sin'
factor_b = 1

# --------------------------- Extraigo datos
# 1.7933
datos = DatosProcesadosDiff(path, factor_b=factor_b,
                            bmax=bmax, gpshape=gpshape)
delta = datos.delta
bigDelta = datos.bigDelta
gpmax = datos.gpmax
# -----------------------------------------------


# -----------------------------------------------
# datos de la muestra
muestra = f"M{Nmuestra}"
N = Nmuestra-10
# matriz es el medio: Bulk, Q30 o Q3
if Nmuestra == 999:
  matriz = 'bulk'
elif N > 10:
  if Nmuestra == 21:
    matriz = 'Q3'
  elif Nmuestra == 22:
      matriz = 'Q30'
  elif Nmuestra == 24:
      matriz = 'Bulk'
elif Nmuestra==0:
    matriz = 'referencia LiCl 0.25 M'
elif N < 3:
    matriz = 'Q3'
elif N < 6:
    matriz = 'Q30'
elif N < 9:
    matriz = 'Bulk'
# pc es el porcentaje de glicerol: 50%, 70% o 90%
if Nmuestra == 999:
  pc = 0
elif N > 10:
    pc = 30
elif Nmuestra==0:
    pc = 0
elif N % 3 == 0:
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
         fr"$\Delta = {bigDelta}$ ms, $\delta = {delta}$ ms, "\
         fr"$G_{{max}} = {gpmax:.2f} \%$"
axs[0].set_title(titulo)
axs[0].plot(bvalue, signal, 'ko')
axs[0].plot(bvalue_fit, signal_fit, 'r-')
text = f"$D =$ {datos.fit_results[0]:.4f} $10^{{-9}}m^2/s$ \n " \
       f"$u_D =$ {datos.fit_results[1]:.1g} $10^{{-9}}m^2/s$ \n " \
       f"$r^2 =$ {datos.fit_results[2]:.5g}"
axs[0].text(bvalue[-1]*0.6, 0.7*signal[0], text,
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
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}/datos_Diff/figuras/{filename0}_Diff.png'
    fig.savefig(filename)   # save the figure to file

    filename = f'{savepath}/datos_Diff/figuras/' \
               f'{filename0}_Diff-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n"\
             f"Rango de integracion: {ppmRange} ppm\n"\
             f"bvalue (10^-9 m^2/s)\t S (norm)"
    Diffdata = np.array([bvalue, signal]).T
    np.savetxt(f"{savepath}/datos_Diff/{filename0}_Diff.dat",
               Diffdata, header=header)

    data = np.array([ppmAxis, re, im]).T
    filename = f"{savepath}/datos_Diff/primerEspectro/" \
               f"{filename0}_primerEspectro.dat"
    np.savetxt(filename, data)
