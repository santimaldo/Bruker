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


# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"

# path_bruker = "2025-06-17_3.2mm_IMECdendrites___cal/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025_06_IMEC/"
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


# savepath_especifico = ""
# nucleo = "1H"
# # info: muestra, expn, ppmRange
# infos = [['hanecho', 100, None],
#          ['cp_mwON_bsms_-2400_p15_1500us', 37, None],
#          ['cp_mwON_bsms_9999_p15_1500us', 55, None],
#          ['cp_mwOFF_bsms_-2400_p15_1500us', 101, None],
#          ]

# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# infos = [#['hanecho', 100, None],
#          ['cp_mwON_bsms_-2400_p15_1500us', 102, None],
#          #['cp_mwON_bsms_9999_p15_1500us', 55, None],
#           ['zg_mwON_bsms_-2400_p15_1500us', 101, None],
#           ['zg_mwOFF_bsms_-2400_p15_1500us', 3, None],
#          ]




path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "2025-06-27_InSitu/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\Li-on-mesh"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [550, -150]  # rango de ppm a integrar
infos = [
         ['zg_mwOFF_bsms_-2000_BYPASS', 1, ppmRange],
         ['zg_mwON_bsms_-2000_190K', 6, ppmRange],
         ['zg_mwON_bsms_-3550_190K', 18, ppmRange],
         ['zg_mwON_bsms_9999_190K', 24, ppmRange],
         ['zg_mwOFF_bsms_-2000_290K', 30, ppmRange],
         ]

path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "2025-06-27_InSitu/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\Li-on-mesh-LFP"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [550, -150]  # rango de ppm a integrar
infos = [
         ['zg_mwOFF_bsms_0_2.4A_290K', 100, ppmRange],
         ['zg_mwOFF_bsms_0_2.4A_190K', 107, ppmRange],
         ['zg_mwON_bsms_1250_2.2A_190K', 135, ppmRange],
         ['zg_mwON_bsms_-9999_2.2A_190K', 141, ppmRange],
         ]

path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "2025-06-26_InSitu/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\LiDendrites-and-KBr"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [550, -150]  # rango de ppm a integrar
infos = [
         ['hahnecho_mwON_bsms_9999_2.2A_190K', 22, ppmRange],
         ['hahnecho_mwON_bsms_-2000_2.4A_190K', 34, ppmRange],
         ['hahnecho_mwOFF_bsms_-2000_2.4A_190K', 38, ppmRange],
         ['hahnecho_mwOFF_bsms_-2000_2.4A_290K', 40, ppmRange],
         ]


path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "InSitu_February_2025/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-perpendicular/"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [300, -500]  # rango de ppm a integrar
infos = [
         ['OFF', 1, ppmRange],
         ['ON', 8, ppmRange],
         ['offNOE', 12, ppmRange]
         ]

path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "InSitu_February_2025/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-LFP-perpendicular/"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [300, -500]  # rango de ppm a integrar
infos = [
         ['OFF', 16, ppmRange],
         ['ON', 20, ppmRange],
         ['offNOE', 24, ppmRange]
         ]

path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "InSitu_February_2025/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-LFP-parallele/"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [300, -500]  # rango de ppm a integrar
infos = [
         ['OFF', 28, ppmRange],
         ['ON', 34, ppmRange],
         ['offNOE', 39, ppmRange]
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
    ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    # ax.set_xlim([50,-50])
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