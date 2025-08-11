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


path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500/"

path_bruker = "2025-05-15_PEO-solid-electrolyte/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2025-05_500MHz_PEO_PEO-PTT/"
savepath_especifico = ""
# info: muestra, expn, ppmRange
nucleo = "1H"
ppmRange = None
infos = [['PEO-LiTFSI_6kHz', 3],
         ['PEO-LiTFSI_10kHz', 4],
         ['PEO-LiTFSI_15kHz', 12],
         ['PEO-LiTFSI_20kHz', 13],
         ['PEO-LiTFSI_25kHz', 14],
         ['PEO-LiTFSI_30kHz', 15],
         
         ['PEO-PTT-LiTFSI_20kHz', 46],
         ['PEO-PTT-LiTFSI_30kHz', 84]
         ]


# info: muestra, expn, ppmRange
nucleo = "7Li"
ppmRange = None
infos = [['PEO-LiTFSI_10kHz', 25], # compare with 7 --> before fast spinning
         ['PEO-LiTFSI_20kHz', 23],
         ['PEO-LiTFSI_30kHz', 19],
         
         ['PEO-PTT-LiTFSI_10kHz', 34],
         ['PEO-PTT-LiTFSI_15kHz', 40],
         ['PEO-PTT-LiTFSI_20kHz', 42],
         ['PEO-PTT-LiTFSI_25kHz', 80],
         ['PEO-PTT-LiTFSI_30kHz', 82]
         ]

# info: muestra, expn, ppmRange
nucleo = "19F"
ppmRange = None #[-85,-75]
infos = [['PEO-LiTFSI_10kHz', 26], # compare with 7 --> before fast spinning
         ['PEO-LiTFSI_20kHz', 24],
         ['PEO-LiTFSI_30kHz', 21],
         
         ['PEO-PTT-LiTFSI_10kHz', 36],
         ['PEO-PTT-LiTFSI_15kHz', 37],
         ['PEO-PTT-LiTFSI_20kHz', 41],
         ['PEO-PTT-LiTFSI_25kHz', 81],
         ['PEO-PTT-LiTFSI_30kHz', 83]
         ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-07-16_PEO-solid-electrolyte_VT/"
# path_bruker = ""
# savepath_local = r""
# savepath_especifico = ""
# nucleo = "7Li"
# ppmRange = None #[-85,-75]
# infos = [['tmp', 111]]


fig, ax = plt.subplots(num=1, nrows=1, ncols=1)  # create figure & 1 axis
for info in infos:
    muestra, expn = info
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