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

path_bruker = "InSitu_May_2025/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-05_mesh/"
# savepath_especifico = "spec1d/"
# # info: muestra, expn, ppmRange
# infos = [['uW-ON_mesh', 31, None],
#          ['uW-ON-off-field_mesh', 33, None],
#          ['uW-OFF--cool_mesh', 23, None],
#          ['uW-OFF--warm_mesh_AFTER', 35, None],
#          ['uW-OFF--warm_mesh_BEFORE', 11, None]
#          ]


savepath_especifico = "spec1d-LFP/"
# info: muestra, expn, ppmRange
infos = [['uW-ON_mesh-lfp', 58, None],
         ['uW-ON-off-field_mesh-lfp', 59, None],
         ['uW-OFF--warm_mesh_AFTER-lfp', 81, None],
         ['uW-OFF--warm_mesh_BEFORE-lfp', 51, None]
         ]

# path_bruker = "3.2mm-Santi-IMECdendrites-2025-03-21/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\analysis\2025-03_IMEC_sample/"
# savepath_especifico = ""
# # info: muestra, expn, ppmRange
# infos = [['uW-OFF_o1-metal', 17, None],
#          ['uW-ON_o1-metal', 9, None],
#          ['uW-ON-offset_o1-metal', 12, None],
#          ['uW-OFF_o1-diamgnetic', 19, None],
#          ['uW-ON_o1-diamagnetic', 11, None],
#          ['uW-ON-offset_o1-diamagnetic', 13, None],
#          ]

nucleo = "7Li"
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