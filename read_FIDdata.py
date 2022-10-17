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
info = [16, '09-29', 1, [-1, 1], 6]
info = [17, '09-29', 8, [-1, 1], np.inf]
# info = [18, '10-03', 10, [-1, 1], 80]
# # Q30:
# info = [13, '09-29', 4, [-1.3, 0.7], 40]
# info = [14, '10-03', 1, [-1.5, 0.8], 30]
# info = [15, '10-11', 1, [-2.4, 1.17], np.inf]
# Q3:
info = [10, '10-11', 12, [-1.3, 0.7], 40]


save = True
Nmuestra, fecha, expn, ppmRange, bmax = info
#-------------------- directorios
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_local = "S:/PosDoc/Glicerol-Agua/116MHz"

path_bruker = f"/2022-{fecha}_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/"
savepath = "S:/temp"

# --------------------------- Extraigo datos
datos = Datos(path, set_fid=True)
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
timeAxis = datos.fid.timeAxis
# spec = datos.espectro.real

re = datos.fid.real
im = datos.fid.imag

titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n"

fig1d, ax1d = plt.subplots(1, 1)
ax1d.axhline(0, color='gray')
ax1d.plot(timeAxis, re/np.max(re), 'k', lw=2)
ax1d.plot(timeAxis, im/np.max(re), 'r', lw=2)
# ax1d.set_xlim(np.max(timeAxis), np.min(timeAxis))
ax1d.set_xlabel("NMR Shift [ppm]")
ax1d.set_title(titulo)

plt.axhline(0, color='k')
# -----------------------------------------------

# %%
# guardo data:
if save:
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}/figuras/{filename0}_FID.png'
    fig1d.savefig(filename)   # save the figure to file


    # header = f"Archivo de datos: {path_bruker}\n"\
    #          f"Rango de integracion: {ppmRange} ppm\n"\
    #          f"bvalue (10^-9 m^2/s)\t S (norm)"
    # Diffdata = np.array([bvalue, signal]).T
    # np.savetxt(f"{savepath}/{filename0}_Diff.dat",
    #            Diffdata, header=header)

    # data = np.array([timeAxis, re, im]).T
    # filename = f"{savepath}/primerEspectro/" \
    #            f"{filename0}_primerEspectro.dat"
    # np.savetxt(filename, data)
