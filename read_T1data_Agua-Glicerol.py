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


# Nmuestra es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
Nmuestra = 15
# directorio de datos:
fecha = '10-11'  # 'MM-DD'
expn = 2
# rango de integracion
ppmRange = [-1.5, 1]
save = True

#-------------------- directorios
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_bruker = f"/2022-{fecha}_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/datos/"


# --------------------------- Extraigo datos
datos = DatosProcesadosT1(path)
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

re = datos.espectro.real[-1]
im = datos.espectro.imag[-1]

titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n"

r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
fig1d, ax1d = plt.subplots(1, 1, num=7531)
ax1d.axvspan(r1, r2, alpha=0.2, color='b')
ax1d.axhline(0, color='gray')
ax1d.plot(ppmAxis, re/np.max(re), 'k', lw=2)
ax1d.text(r1-1, 0.8, "Region de integracion (T1)", color='b')
ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax1d.set_xlabel("NMR Shift [ppm]")

ax1d.set_title(titulo)

plt.axhline(0, color='k')
# -----------------------------------------------
# %%

tau, signal = datos.get_T1data(ppmRange)
tau_fit, signal_fit, residuals = datos.T1fit()

titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n"

fig, axs = plt.subplots(2, 2, figsize=(8, 6))
fig.suptitle(titulo)
# -------------
axs[0, 0].plot(tau, signal, 'ko')
axs[0, 0].plot(tau_fit, signal_fit, 'r-')
text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
               multialignment="left")
axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
# -------------
axs[1, 0].plot(tau, residuals, 'ko')
axs[1, 0].axhline(0, color='k', linestyle='--')
axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
# -------------
axs[0, 1].plot(tau, signal, 'ko')
axs[0, 1].plot(tau_fit, signal_fit, 'r-')
axs[0, 1].set(xlabel=r'$\tau$ [ms]')
axs[0, 1].set_xscale('log')
# -------------
axs[1, 1].plot(tau, residuals, 'ko')
axs[1, 1].axhline(0, color='k', linestyle='--')
axs[1, 1].set(xlabel=r'$\tau$ [ms]')
axs[1, 1].set_xscale('log')


for ax in axs.flat:
    ax.label_outer()


# guardo data:
if save:
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}/{filename0}_T1.png'
    fig.savefig(filename)   # save the figure to file

    filename = f'{savepath}/{filename0}_T1-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n"\
             f"Rango de integracion: {ppmRange} ppm\n"\
             f"tau (ms)\t S (norm)"
    T1data = np.array([tau, signal]).T
    np.savetxt(f"{savepath}/{filename0}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{filename0}_T1-ultimoEspectro.dat", data)
