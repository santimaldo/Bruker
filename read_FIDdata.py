# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Espectro import autophase
import scipy.integrate as integrate
import matplotlib.ticker as ticker


# Nmuestra:  es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
# fecha: correspondiente al directorio de datos.  'MM-DD'
# expn: numero de archivo
# info = [Nmuestra, fecha, expn]
# Bulk:
info = [16, '09-29', 1]
info = [17, '09-29', 8]
info = [18, '10-03', 10]
# # # Q30:
info = [13, '09-29', 4]
info = [14, '10-03', 1]
info = [15, '10-11', 1]
# Q3:
info = [10, '10-03', 24]


save = False
Nmuestra, fecha, expn = info
#-------------------- directorios
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
# path_local = "S:/PosDoc/Glicerol-Agua/116MHz"

path_bruker = f"/2022-{fecha}_Diff_Silica_Agua-Glicerol-LiCl/{expn}/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/7Li/datos_FID/"


# --------------------------- Extraigo datos
datos = Datos(path, set_fid=True)
timeAxis = datos.fid.timeAxis
# spec = datos.espectro.real

re = datos.fid.real
im = datos.fid.imag
mod = np.abs(re+1j*im)

fid, angle = autophase(re+1j*im, x=timeAxis, method='maxIntReal_fid')

re = fid.real
im = fid.imag
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


titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n"

fig1d, ax1d = plt.subplots(1, 1)
ax1d.axhline(0, color='gray')
ax1d.plot(timeAxis*1000, re, 'k', lw=2, label='Real')
ax1d.plot(timeAxis*1000, im, 'r', lw=2, label='Imag')
ax1d.plot(timeAxis*1000, mod, 'orange', lw=2, alpha=0.4, label='Abs')
# ax1d.set_xlim(np.max(timeAxis), np.min(timeAxis))
# ax1d.plot(timeAxis, np.real(fid), 'k--', lw=2, label='Real')
# ax1d.plot(timeAxis, np.imag(fid), 'r--', lw=2, label='Imag')

ax1d.set_xlabel("Time [ms]")
ax1d.set_title(titulo)
ax1d.legend()
plt.axhline(0, color='k')
# -----------------------------------------------

# %%
# guardo data:
if save:
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}figuras/{filename0}_FID.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n"\
             f"time (s)\t Real (raw)\t Imag (raw)\t Abs (raw)"
    Diffdata = np.array([timeAxis, re, im, mod]).T
    np.savetxt(f"{savepath}/{filename0}_FID.dat",
               Diffdata, header=header)
