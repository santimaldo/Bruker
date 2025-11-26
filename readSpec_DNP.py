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

# path_bruker = "2025-06-17_3.2mm_IMECdendrites___cal/"
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


# savepath_especifico = ""
# nucleo = "1H"
# # info: muestra, expn, ppmRange
# infos = [['hanecho', 100, None],
#          ['cp_mwON_bsms_-2400_p15_1500us', 37, None],
#          ['cp_mwON_bsms_9999_p15_1500us', 55, None],
#          ['cp_mwOFF_bsms_-2400_p15_1500us', 101, None],
#          ]


path_bruker = "2025-06-17_3.2mm_IMECdendrites___cal/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DegradationProject\2025-07-consortium/"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
infos = [#['hanecho', 100, None],
         ['eLi_expn102_cp_mwOFF_bsms_-2400_p15_1500us', 102, None],
         #['cp_mwON_bsms_9999_p15_1500us', 55, None],
          ['eLi_expn31_zg_mwON_bsms_-2400', 31, None],
          ['eLi_expn3_zg_mwOFF_bsms_-2400', 3, None],
         ]

path_bruker = "3.2mm-Santi-IMECdendrites-2025-04-28___cal/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DegradationProject\2025-07-consortium/"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
infos = [#['hanecho', 100, None],
          ['LiOH_expn40_zg_mwOFF', 40, None],
          ['LiOH.H20_expn30_zg_mwOFF', 30, None],
         ]

path_bruker = "2025-06-17_3.2mm_IMECdendrites___cal/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DegradationProject\2025-07-consortium/1H/"
savepath_especifico = ""
nucleo = "1H"
# info: muestra, expn, ppmRange
infos = [#['hanecho', 100, None],
         ['eLi_cp_mwOFF_bsms_-2400_p15_1500us', 101, None],
         ['eLi_cp_mwON_bsms_-2400_p15_1500us', 37, None],
         ['eLi_hahnecho_mwOFF_bsms_-2400', 100, None],
         ['eLi_cp_mwON_bsms_-2400_p15_2000us', 2199903, None],
         ['eLi_cp_mwON_bsms_-2400_p15_150us', 2199914, None],
         ]

# path_bruker = "3.2mm-Santi-IMECdendrites-2025-04-28___cal/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DegradationProject\2025-07-consortium/"
# savepath_especifico = ""
# nucleo = "1H"
# # info: muestra, expn, ppmRange
# infos = [#['hanecho', 100, None],
#           ['LiOH_hahnecho_mwOFF', 41999, None],
#           ['LiOH.H20_hahnecho_mwOFF', 31, None],
#          ]


# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-06-27_InSitu/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\Li-on-mesh"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [550, -150]  # rango de ppm a integrar
# infos = [
#          ['zg_mwOFF_bsms_-2000_BYPASS', 1, ppmRange],
#          ['zg_mwON_bsms_-2000_190K', 6, ppmRange],
#          ['zg_mwON_bsms_-3550_190K', 18, ppmRange],
#          ['zg_mwON_bsms_9999_190K', 24, ppmRange],
#          ['zg_mwOFF_bsms_-2000_290K', 30, ppmRange],
#          ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-06-27_InSitu/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\Li-on-mesh-LFP"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [550, -150]  # rango de ppm a integrar
# infos = [
#          ['zg_mwOFF_bsms_0_2.4A_290K', 100, ppmRange],
#          ['zg_mwOFF_bsms_0_2.4A_190K', 107, ppmRange],
#          ['zg_mwON_bsms_1250_2.2A_190K', 135, ppmRange],
#          ['zg_mwON_bsms_-9999_2.2A_190K', 141, ppmRange],
#          ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-06-26_InSitu/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\LiDendrites-and-KBr"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [550, -150]  # rango de ppm a integrar
# infos = [
#          ['hahnecho_mwON_bsms_9999_2.2A_190K', 22, ppmRange],
#          ['hahnecho_mwON_bsms_-2000_2.4A_190K', 34, ppmRange],
#          ['hahnecho_mwOFF_bsms_-2000_2.4A_190K', 38, ppmRange],
#          ['hahnecho_mwOFF_bsms_-2000_2.4A_290K', 40, ppmRange],
#          ]


# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "InSitu_February_2025/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-perpendicular/"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [300, -500]  # rango de ppm a integrar
# infos = [
#          ['OFF', 1, ppmRange],
#          ['ON', 8, ppmRange],
#          ['offNOE', 12, ppmRange]
#          ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "InSitu_February_2025/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-LFP-perpendicular/"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [300, -500]  # rango de ppm a integrar
# infos = [
#          ['OFF', 16, ppmRange],
#          ['ON', 20, ppmRange],
#          ['offNOE', 24, ppmRange]
#          ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "InSitu_February_2025/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_05_Quartz\February_Marie_Mesh-LFP-parallele/"
# savepath_especifico = ""
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# ppmRange = [300, -500]  # rango de ppm a integrar
# infos = [
#          ['OFF', 28, ppmRange],
#          ['ON', 34, ppmRange],
#          ['offNOE', 39, ppmRange]
#          ]

########## 31 julio 2025 - Li on Cu mesh ##########
path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "2025-07-31_InSitu/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_08_Quartz\Li-on-mesh"
savepath_especifico = ""
nucleo = "7Li"
# info: muestra, expn, ppmRange
ppmRange = [550, -150]  # rango de ppm a integrar
infos = [
         ['zg_mwON_bsms_3850_2.2A_190K', 17, ppmRange],
         ['zg_mwON_bsms_3100_2.2A_190K', 18, ppmRange],
         ['zg_mwON_bsms_-9999_2.2A_190K', 22, ppmRange],
         ['zg_mwOFF_bsms_0_2.2A_190K', 1, ppmRange],
         ]
# ### 31 julio 2025 - Li on Cu mesh - Al mesh ######
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-07-31_InSitu/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025_08_Quartz\Li-on-mesh-Al-mesh/"
# savepath_especifico = ""
# nucleo = "7Li"
# info: muestra, expn, ppmRange
# ppmRange = [550, -150]  # rango de ppm a integrar
# infos = [
#          ['zg_mwOFF_bsms_0_2.2A_RT', 30, ppmRange],
#          ['zg_mwON_bsms_0_2.2A_190K', 33, ppmRange],
#          ['zg_mwON_bsms_-6000_2.4A_190K', 54, ppmRange],
#          ]




#### 27 Oct 2025 - Li on CU Mesh + PTFE holders --- Bruker
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-10-27_InSitu/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-10_InSitu\LionCuMesh/"

# ppmRange = [1000, -200]  # rango de ppm a integrar
# nucleo = "7Li"
# savepath_especifico = "Angle_000deg"
# # info: muestra, expn, ppmRange
# infos = [
#          ['mwOFF_2.2A_bsms_0', 2, ppmRange],
#          ['mwON_2.2A_bsms_9999_ultra-off-OE', 15, ppmRange],
#          ['mwON_2.5A_bsms_9999_off-OE', 26, ppmRange],
#          ['mwON_2.2A_bsms_-4050_OE-metal', 22, ppmRange],
#          ['mwON_2.2A_bsms_-4450_OE-SEI', 23, ppmRange],
#          ]
# savepath_especifico = "Angle_180deg"
# # info: muestra, expn, ppmRange
# infos = [
#          ['mwOFF_2.5A_bsms_0', 40, ppmRange],
#          ['mwON_2.5A_bsms_50_OE-metal', 48, ppmRange],
#          ['mwON_2.5A_bsms_-200_OE-SEI', 49, ppmRange],
#          ]
# savepath_especifico = "Angle_090deg"
# # info: muestra, expn, ppmRange
# infos = [
#          ['mwOFF_2.5A_bsms_-200', 60, ppmRange],
#          ['mwON_2.5A_bsms_50_OE-metal', 70, ppmRange],
#          ['mwON_2.5A_bsms_-200_OE-SEI', 69, ppmRange],
#          ['mwOFF_2.5A_bsms_-200_next-day', 71, ppmRange]
#          ]


# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-10_InSitu\LionCuMesh_LFP/"
# ppmRange = [800, -200]  # rango de ppm a integrar
# nucleo = "7Li"
# savepath_especifico = ""
# # info: muestra, expn, ppmRange
# infos = [
#          ['mwOFF_2.5A_bsms_0', 100, ppmRange],
#          ['mwON_2.5A_bsms_9999_off-OE', 112, ppmRange],
#          ['mwON_2.2A_bsms_-4050_OE-metal', 111, ppmRange]
#          ]








# ### 30 sep 2025 - Li polymer + KBr ---Kieran
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-09-30_Kieran_Li-polymer-KBr_cal/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\Kieran/"
# savepath_especifico = ""
# ppmRange = [500, -150]  # rango de ppm a integrar
# nucleo = "7Li"
# # info: muestra, expn, ppmRange
# infos = [
#          ['7Li_mwON_bsms_-3500_2.2A', 20, ppmRange],
#          ['7Li_mwON_bsms_-3500_2.2A_o1diamagnetic', 21, ppmRange],
#          ['7Li_mwOFF_bsms_-3500_2.2A', 30, ppmRange]
#          ]
# nucleo = "6Li"
# # info: muestra, expn, ppmRange
# infos = [
#          ['6Li_mwON_bsms_-3500_2.2A', 24, ppmRange],
#          ['6Li_mwOFF_bsms_-3500_2.2A', 31, ppmRange]
#          ]

# ####### Kieran 2025-10
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
# path_bruker = "2025-10-24_Kieran_Li-polymer-KBr_cal/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\Kieran\2025-09_without-polymer/"

# ppmRange = [800, -200]  # rango de ppm a integrar
# # nucleo = "7Li"
# # savepath_especifico = ""
# # # info: muestra, expn, ppmRange
# # infos = [
# #          ['mwOFF_2.2A_bsms_0', 1, ppmRange],
# #          ['mwON_2.2A_bsms_-6500_o1metal', 8, ppmRange],
# #          ['mwON_2.2A_bsms_-6500_o1diamagnetic', 14, ppmRange]
# #          ]
# nucleo = "6Li"
# # info: muestra, expn, ppmRange
# infos = [
#          ['6Li_mwON_bsms_-6500_2.2A', 9, ppmRange]
#         #  ['6Li_mwOFF_bsms_-6500_2.2A', 31, ppmRange]
#          ]

#=============================================================
####### Rui 2025-11
path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp/"
path_bruker = "2025-11-13_3.2mm_Rui-dendrites/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-11_DNP_CC\spec/"

nucleo = "7Li"
savepath_especifico = ""
## info: muestra, expn, ppmRange
infos = [
         ['mwOFF_2.2A_bsms_-1600_o1metal', 105, [400, 100]],
         ['mwON_2.2A_bsms_-1600_o1metal', 83, [400,100]],
         ['mwOFF_2.2A_bsms_-1600_o1diamagnetic', 103, [150,-150]],
         ['mwON_2.2A_bsms_-1600_o1diamagnetic', 79, [150,-150]],
         ['mwOFF_2.2A_bsms_-1600_LiF-reference', 110, [150,-150]],
         ['mwON-offOE_2.2A_bsms_9999_o1metal', 58, [400, 100]],
         ['mwON-offOE_2.2A_bsms_9999_o1diamagnetic', 59, [150,-150]]
         ]
# ppmRange = [100, -500]  # rango de ppm a integrar
# nucleo = "19F"
# savepath_especifico = ""
# ## info: muestra, expn, ppmRange
# infos = [
#          ['mwOFF_2.2A_bsms_-1600_CP', 102, ppmRange],
#          ['mwON_2.2A_bsms_-1600_CP', 76999, ppmRange],
#          ['mwON-offOE_2.2A_bsms_9999_CP', 77, ppmRange],
#          ['mwOFF_2.2A_bsms_-1600_LiF-reference', 111, ppmRange]         
#          ]

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