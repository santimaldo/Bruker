# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:09:33 2022

@author: Santi

Extrae multiples espectros adquiridos en Bruker Avance II
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Espectro import autophase

path = 'S:/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/'
savepath = 'S:/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/Datos/'

expn = 104
muestra = "LiTf-G2_M4_HF_Topspin_lb10"
nucleo = "7Li"

save = True

ppmRange = [-30,30]


# extraigo:
datos = DatosProcesados(f'{path}{expn}/')
datos.espectro.ppmSelect(ppmRange)
re = datos.espectro.real
im = datos.espectro.imag
re_norm = re/np.max(re)
im_norm = im/np.max(re)
ppmAxis = datos.espectro.ppmAxis



# grafico para guardar:
fig, ax = plt.subplots(num=1, nrows=1, ncols=1 )  # create figure & 1 axis
ax.plot(ppmAxis, re_norm, linewidth=2)
ax.axhline(0, color='k', ls='--')
ax.set_title(muestra)
ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])

if save:
  # guardo:
  header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
  dataexport = np.array([ppmAxis, re_norm, im_norm, re, im]).T
  filename = f'{savepath}/{nucleo}_{muestra}.dat'
  np.savetxt(filename, dataexport, header=header)

  filename = f'{savepath}/{nucleo}_{muestra}.png'
  fig.savefig(filename)   # save the figure to file
  plt.close(fig)    # close the figure window

plt.show()
