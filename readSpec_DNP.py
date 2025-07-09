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
from VoigtFit import VoigtFit


path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"

path_bruker = "2025-06-17_3.2mm_IMECdendrites/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025_06_IMEC/"


# savepath_especifico = "CP_7Li-1H_uW-OFF/"
# nucleo = "1H"
# # info: muestra, expn, ppmRange
# infos = [['uwOFF_p15-3000us', 110, None],
#          ['uwOFF_p15-1000us', 111, None],
#          ['uwOFF_p15-0500us', 112, None],
#          ['uwOFF_p15-0400us', 113, None],
#          ['uwOFF_p15-0300us', 114, None],
#          ['uwOFF_p15-0200us', 115, None],
#          ['uwOFF_p15-0100us', 116, None],
#          ['uwOFF_p15-0050us', 117, None],
#          ]


savepath_especifico = ""
nucleo = "1H"
# info: muestra, expn, ppmRange
infos = [['hanecho', 100, None],
         ['cp_mwON_bsms_-2400_p15_1500us', 37, None],
         ['cp_mwON_bsms_9999_p15_1500us', 55, None],
         ]


fig, ax = plt.subplots(num=1, nrows=1, ncols=1)  # create figure & 1 axis
for info in infos:
    muestra, expn, ppmRange = info
    save = True


    # extraigo:
    datos = DatosProcesados(f'{path_local}{path_bruker}{expn}/')
    # datos = DatosProcesados(f'{path_local}')
    if ppmRange is not None:
        datos.espectro.ppmSelect(ppmRange)
    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis


    # grafico para guardar:
    ax.plot(ppmAxis, re, linewidth=2)
    ax.axhline(0, color='k', ls='--')
    ax.set_title(muestra)
    ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
    # ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    ax.set_xlim([50,-50])
    if save:
        savepath = f"{savepath_local}{savepath_especifico}"
        # guardo:
        header = "ppmAxis\t real \t imag"
        dataexport = np.array([ppmAxis, re, im]).T
        filename = f'{savepath}/{nucleo}_Nexp{expn}_{muestra}.dat'
        np.savetxt(filename, dataexport, header=header)

        # filename = f'{savepath}/{nucleo}_{muestra}.png'
        # fig.savefig(filename)   # save the figure to file
        # # plt.close(fig)    # close the figure window

plt.show()


# vfit=VoigtFit(ppmAxis, 
#               re, 
#               Npicos=4,
#               ajustar=True,
#               center=[0, 0,-1,-3]
#               )
# fig = vfit.plot_ajuste()
# fig.gca().set_xlim([8,-8])