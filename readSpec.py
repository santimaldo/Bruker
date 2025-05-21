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


path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/"
path_bruker = "400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP/"
savepath_especifico = "1Dspec/bsms_field_-5850/"



sample = "IMEC"

# info: nucleus,muestra, expn, ppmRange
# info = ["1H", "IMECdendrites_KBr", 52, [350,100]]
# info = ["1H", "LiOH", 48, [350,100]]
# info = ["1H", "LiOH.H20", 31, [350,100]]
# info = ["1H", "LiOH_after-one-day", 41999, [350,100]]
#info = ["1H", "eLi_CP_ct-01ms_d1-2s", 1699801, [350,100]]
# info = ["1H", "eLi_CP_ct-10ms_d1-2s", 1699601, [350,100]]
# info = ["1H", "eLi_uW-OFF", 26, [350,100]]
# info = ["1H", "eLi_uW-ON", 17, [350,100]]
# reference_value = -1.5 # ppm    
# reference_measured = 226.7633
########### 7Li
# info = ["7Li", "eLi_uW-OFF", 55, [-225, -320]]
# info = ["7Li", "eLi_uW-ON", 19, [-225, -320]]
# info = ["7Li", "LiOH", 40, [-225, -320]]
# info = ["7Li", "LiOH.H20", 30, [-225, -320]]
# reference_value = 0.4# ppm    
# reference_measured = -273.9571081
# - - - - - - no baselinsubstraction
# info = ["7Li", "eLi_uW-OFF", 55, [-114, -600]]
# info = ["7Li", "eLi_uW-ON", 19, [-114, -600]]
# # info = ["7Li", "LiOH", 40, [100, -600]]
# # info = ["7Li", "LiOH.H20_no-basline-subs", 30, [100, -600]]
info = ["7Li", "eLi_o1-metal_uW-OFF", 25, [200,-500]]
# info = ["7Li", "eLi_o1-metal_uW-ON", 13, [200,-500]]
reference_value = 0.4# ppm    
reference_measured = -273.9

ppmCorrection = reference_value - reference_measured
nucleo, muestra, expn, ppmRange = info
save = True


# extraigo:
datos = DatosProcesados(f'{path_local}{path_bruker}{expn}/')
# datos = DatosProcesados(f'{path_local}')
if ppmRange is not None:
    datos.espectro.ppmSelect(ppmRange)
re = datos.espectro.real
im = datos.espectro.imag
re_norm = re/np.max(re)
im_norm = im/np.max(re)
ppmAxis = datos.espectro.ppmAxis
ppmAxis = ppmAxis + ppmCorrection  # ppm axis correction

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


# vfit=VoigtFit(ppmAxis, 
#               re, 
#               Npicos=4,
#               ajustar=True,
#               center=[0, 0,-1,-3]
#               )
# fig = vfit.plot_ajuste()
# fig.gca().set_xlim([8,-8])