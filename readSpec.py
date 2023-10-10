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

path_local = "S:/NMRdata/2023_InSitu/"
path_bruker = "2023-05-12_LiMetal_in-situ/"

savepath_local = "S:/Posdoc/InSitu/"
savepath_especifico = "2023-05_test1/"


path_local = "S:/Doctorado/LiMetal/116MHz/"
path_bruker = "2022-03-23_SMC_test/54/"
savepath_local = ""
savepath_especifico = "S:/tmp/"

# info: muestra, expn, ppmRange
info = ['t11_ciclado_10mA', 8, None]

nucleo = "7Li"
muestra, expn, ppmRange = info

expn = ""
save = True


# extraigo:
datos = DatosProcesados(f'{path_local}{path_bruker}{expn}/')
if ppmRange is not None:
    datos.espectro.ppmSelect(ppmRange)
re = datos.espectro.real
im = datos.espectro.imag
re_norm = re/np.max(re)
im_norm = im/np.max(re)
ppmAxis = datos.espectro.ppmAxis


# grafico para guardar:
fig, ax = plt.subplots(num=1, nrows=1, ncols=1)  # create figure & 1 axis
ax.plot(ppmAxis, re_norm, linewidth=2)
ax.axhline(0, color='k', ls='--')
ax.set_title(muestra)
ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])

if save:
    savepath = f"{savepath_local}{savepath_especifico}"
    # guardo:
    header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
    dataexport = np.array([ppmAxis, re_norm, im_norm, re, im]).T
    filename = f'{savepath}/{nucleo}_{muestra}.dat'
    np.savetxt(filename, dataexport, header=header)

    filename = f'{savepath}/{nucleo}_{muestra}.png'
    fig.savefig(filename)   # save the figure to file
    # plt.close(fig)    # close the figure window

plt.show()
