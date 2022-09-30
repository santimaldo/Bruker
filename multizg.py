# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:09:33 2022

@author: Santi

Extrae multiples espectros adquiridos en Bruker Avance II
"""

# import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

#----------------------------------------------------
path = 'S:/Carbones/300MHz/2022-09-01_CarbonesCNEA_M4_LiTFSI/'
savepath = 'S:/Carbones/analisis/2022-09_M4_LiTFSI-G2_100mM/files/'

expnums = [[6, 100, 200],     #  1H
           [7, 104, 204] ]     #  7Li
muestras = ["LiTFSI-G2_100mM_bulk", "M4-NaOH", "M4-HF"]
nucleos = ["1H", "7Li"]


#-------------------------------------------------------
path = 'S:/Carbones/300MHz/2018-09-25_Carbones_Glyma/'
savepath = 'S:/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/Datos/'

expnums = [#[6, 100, 200],     #  1H
           [16,18] ]     #  7Li
muestras = ["LiTf-G2_100mM_bulk_SHIM_HRMAS", "LiCl_solid_SHIM_HRMAS"]
nucleos = ["7Li"]






ppmRange = [-30,30]

for mm in range(len(nucleos)):
  nucleo = nucleos[mm]
  for nn in range(len(expnums[mm])):
    expn = expnums[mm][nn]
    muestra = muestras[nn]
    print(expn)
    # extraigo:
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(ppmRange)
    re = datos.espectro.real
    im = datos.espectro.imag
    re_norm = re/np.max(re)
    im_norm = im/np.max(re)
    ppmAxis = datos.espectro.ppmAxis

    # guardo:
    header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
    dataexport = np.array([ppmAxis, re_norm, im_norm, re, im]).T
    filename = f'{savepath}/{nucleo}_{muestra}.dat'
    np.savetxt(filename, dataexport, header=header)

    # grafico para ver:
    print('graficando...', nucleo, muestra)
    plt.figure(58)
    plt.plot(ppmAxis, re_norm)
    #plt.plot(ppmAxis, im,'r')
    #plt.plot(mod, 'g')

    # grafico para guardar:
    fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
    ax.plot(ppmAxis, re_norm, linewidth=2)
    ax.set_title(muestra)
    ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
    ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    filename = f'{savepath}/{nucleo}_{muestra}.png'
    fig.savefig(filename)   # save the figure to file
    plt.close(fig)    # close the figure window

plt.show()
