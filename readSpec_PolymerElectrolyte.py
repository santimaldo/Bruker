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


path_local = r"C:\Users\Santi\Documents\NMRdata\500/"

# path_bruker = "2025-05-15_PEO-solid-electrolyte/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2025-05_500MHz_PEO_PEO-PTT/"
# savepath_especifico = ""
# # info: muestra, expn, ppmRange
# nucleo = "1H"
# ppmRange = None
# infos = [['PEO-LiTFSI_6kHz', 3],
#          ['PEO-LiTFSI_10kHz', 4],
#          ['PEO-LiTFSI_15kHz', 12],
#          ['PEO-LiTFSI_20kHz', 13],
#          ['PEO-LiTFSI_25kHz', 14],
#          ['PEO-LiTFSI_30kHz', 15],
         
#          ['PEO-PTT-LiTFSI_20kHz', 46],
#          ['PEO-PTT-LiTFSI_30kHz', 84]
#          ]


# # info: muestra, expn, ppmRange
# nucleo = "7Li"
# ppmRange = None
# infos = [['PEO-LiTFSI_10kHz', 25], # compare with 7 --> before fast spinning
#          ['PEO-LiTFSI_20kHz', 23],
#          ['PEO-LiTFSI_30kHz', 19],
         
#          ['PEO-PTT-LiTFSI_10kHz', 34],
#          ['PEO-PTT-LiTFSI_15kHz', 40],
#          ['PEO-PTT-LiTFSI_20kHz', 42],
#          ['PEO-PTT-LiTFSI_25kHz', 80],
#          ['PEO-PTT-LiTFSI_30kHz', 82]
#          ]

# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = None #[-85,-75]
# infos = [['PEO-LiTFSI_10kHz', 26], # compare with 7 --> before fast spinning
#          ['PEO-LiTFSI_20kHz', 24],
#          ['PEO-LiTFSI_30kHz', 21],
         
#          ['PEO-PTT-LiTFSI_10kHz', 36],
#          ['PEO-PTT-LiTFSI_15kHz', 37],
#          ['PEO-PTT-LiTFSI_20kHz', 41],
#          ['PEO-PTT-LiTFSI_25kHz', 81],
#          ['PEO-PTT-LiTFSI_30kHz', 83]
#          ]

# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-07-16_PEO-solid-electrolyte_VT/"
# path_bruker = ""
# savepath_local = r""
# savepath_especifico = ""
# nucleo = "7Li"
# ppmRange = None #[-85,-75]
# infos = [['tmp', 111]]

#####################################################################################
################# 2026-02 13C CP
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-02_500MHz_13C-CP/"

# path_bruker = "2026-02-06_PEO-solid-electrolyte/"
# savepath_especifico = "PEO-LiTFSI/13C_CP/"
# # info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [100,40]
# infos = [['CP_contactTime_100us', 51, ppmRange],
#          ['CP_contactTime_1000us', 52, ppmRange],
#          ['CP_contactTime_50us', 53, ppmRange],
#          ['CP_contactTime_500us', 54, ppmRange],
#          ['CP_contactTime_5000us', 55, ppmRange],
#          ['CP_contactTime_10us', 56, ppmRange],
#          ['CP_contact-Time_50us_shortD1-0.6s', 60, ppmRange],
#          ['zg_d1-60s', 61, [200, 0]],
#          ]

# path_bruker = "2026-02-07_PEO-PTT_solid-electrolyte/"
# savepath_especifico = "PEO-PTT-LiTFSI/13C_CP/"
# # info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [100,40]
# infos = [['CP_contactTime_100us', 51, ppmRange],
#          ['CP_contactTime_1000us', 52, ppmRange],
#          ['CP_contactTime_50us', 53, ppmRange],
#          ['CP_contactTime_500us', 54, ppmRange],
#          ['CP_contactTime_5000us', 55, ppmRange],
#          ['CP_contactTime_10us', 56, ppmRange],
#          ['CP_contactTime_3000us', 57, ppmRange],
#          ['CP_contact-Time_50us_shortD1-0.6s', 60, ppmRange],
#          ['zg_d1-60s', 61, [200, 0]],
#          ]

# path_bruker = "2026-02-06_PEO-solid-electrolyte/"
# savepath_especifico = "PEO-LiTFSI/7Li/"
# # info: muestra, expn, ppmRange
# nucleo = "7Li"
# ppmRange = [200,-200]
# infos = [['CP_contactTime_2000us_D1-10s', 79, ppmRange],
#          ['CP_contactTime_2000us_D1-0.6s', 78, ppmRange],
#          ]

# path_bruker = "2026-02-06_PEO-solid-electrolyte/"
# savepath_especifico = "PEO-LiTFSI/19F/"
# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = [100,-250]
# infos = [['zg_14kHzMAS', 100, ppmRange],
#          ]

# path_bruker = "2026-02-07_PEO-PTT_solid-electrolyte/"
# savepath_especifico = "PEO-PTT-LiTFSI/19F/"
# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = [100,-250]
# infos = [['PEO-PTT_zg_14kHzMAS', 101, ppmRange],
#          ['PEO-PTT_zg_10kHzMAS', 100, ppmRange]
#          ]

####################################################################################
################ 2026-04 13C CP
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500/"
# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-04_500MHz_13C-CP/"

# path_bruker = "2026-04-19_PEO-PTT_solid-electrolyte/"
# savepath_especifico = "PEO-PTT_LiTFSI/13C_CP/"
# # info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [150,0]
# infos = [['CP_contactTime_500us', 21, ppmRange],
#          ['CP_contactTime_50us', 22, ppmRange],
#          ['CP_contactTime_10us', 23, ppmRange],
#          ['CP_contactTime_100us', 24, ppmRange],
#          ['CP_contactTime_1000us', 25, ppmRange],
#          ['CP_contactTime_3000us', 26, ppmRange],
#          ['CP_contactTime_5000us', 27, ppmRange],
#          ['CP_contactTime_500us_bad-decoupling', 21999, ppmRange],
#          ['CP_contactTime_500us_bad-contact', 21998, ppmRange],
#          ['PEO-PTT_zg-decoupling_d1-60s', 20, [150, 0]],
#          ]

# savepath_especifico = "PEO-PTT/"
# # info: muestra, expn, ppmRange
# nucleo = "1H"
# ppmRange = [100,-100]
# infos = [['PEO-PTT_Hahnecho', 11 , ppmRange]]
# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = [100,-300]
# infos = [['PEO-PTT_zg', 1 , ppmRange]]
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# path_bruker = "2026-04-18_PEO-solid-electrolyte/"
# savepath_especifico = "PEO_LiTFSI/13C_CP/"
# info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [150,0]
# infos = [['CP_contactTime_500us', 21, ppmRange],
#          ['CP_contactTime_50us', 22, ppmRange],
#          ['CP_contactTime_3000us', 23, ppmRange],
#          ['CP_contactTime_100us', 24, ppmRange],
#          ['CP_contactTime_1000us', 25, ppmRange],
#          ['CP_contactTime_10us', 26, ppmRange],
#          ['CP_contactTime_5000us', 27, ppmRange],
#          ['CP_contactTime_500us_bad-parameters', 21999, ppmRange],
#          ['PEO_zg-decoupling_d1-60s', 20, [150, 0]],
#          ]
### - - -  - - -  - - -  - - -  - - -  - - -  - - -  - - -  - - - 
# savepath_especifico = "PEO_LiTFSI/vsT/"
# ####info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [150,0]
# infos = [['CP-after-cooling-back_contactTime_500us', 36, ppmRange],
#          ['CP-after-cooling-back_contactTime_50us', 38, ppmRange],
#         ]
# nucleo = "1H"
# ppmRange = [100,-100]
# infos = [['Hahnecho_antes', 30 , ppmRange],
#          ['Hahnecho_despes', 37 , ppmRange]]
# nucleo = "19F"
# ppmRange = [100,-300]
# infos = [['zg_20deg', 31 , ppmRange],
#          ['zg_30deg', 32 , ppmRange],
#          ['zg_40deg', 33 , ppmRange],
#          ['zg_20deg-after', 39 , ppmRange],]
### - - -  - - -  - - -  - - -  - - -  - - -  - - -  - - -  - - - 
# savepath_especifico = "PEO-PTT/"
# # info: muestra, expn, ppmRange
# nucleo = "1H"
# ppmRange = [100,-100]
# infos = [['PEO-PTT_Hahnecho', 11 , ppmRange]]
# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = [100,-300]
# infos = [['PEO-PTT_zg', 1 , ppmRange]]



####################################################################################
################ 2026-05 VT all nuclei - PEO
savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-05_500MHz_PEO-LiTFSI_VT-phase-trans/"

path_bruker = "2026-05-17_PEO-solid-electrolyte/"
# savepath_especifico = "13C_CP/"
# ## info: muestra, expn, ppmRange
# nucleo = "13C"
# ppmRange = [150,0]
# infos = [['CP_contactTime_100us', 150, ppmRange],
#          ['CP_contactTime_1000us', 151, ppmRange],
#          ['CP_contactTime_10us', 152, ppmRange],
#          ['CP_contactTime_200us', 153, ppmRange],
#          ['CP_contactTime_3000us', 154, ppmRange],
#          ['CP_contactTime_300us', 155, ppmRange],
#          ['CP_contactTime_500us', 156, ppmRange],
#          ['CP_contactTime_50us', 157, ppmRange],
#          ]



savepath_especifico = "13C_CP_vt/"
## info: muestra, expn, ppmRange
nucleo = "13C"
ppmRange = [150,0]
infos = []
for expn in range(32, 133, 10):
    infos.append(["", expn, ppmRange])
infos.append(["", 77, ppmRange])
infos.append(["", 87, ppmRange])
infos.append(["", 97, ppmRange])
infos.append(["", 108, ppmRange])
infos.append(["", 116, ppmRange])

# savepath_especifico = "1H/"
# # info: muestra, expn, ppmRange
# nucleo = "1H"
# ppmRange = [100,-100]
# infos = []
# Temps = [20, 30, 40, 30, 20, 10, 0, -10, 0, 10, 20]
# for expn, temp in zip(range(31, 132, 10), Temps):
#     infos.append([f"Hahnecho_{temp+273}K", expn, ppmRange])
# infos.append([f"Hanecho_{30+273}K_bis", 46, ppmRange])
# infos.append([f"Hanecho_{30+273}K_bis", 66, ppmRange])
# infos.append([f"Hanecho_{20+273}K_bis", 76, ppmRange])
# infos.append([f"Hanecho_{10+273}K_bis", 86, ppmRange])
# infos.append([f"Hanecho_{0+273}K_bis", 96, ppmRange])
# infos.append([f"Hanecho_{-10+273}K_d1-20s", 106, ppmRange])
# infos.append([f"Hanecho_{-10+273}K_d1-30s", 107, ppmRange])

# savepath_especifico = "19F/"
# # info: muestra, expn, ppmRange
# nucleo = "19F"
# ppmRange = [100, -300]
# infos = []
# Temps = [20, 30, 40, 30, 20, 10, 0, -10, 0, 10, 20]
# for expn, temp in zip(range(30, 131, 10), Temps):
#     infos.append([f"Hahnecho_{temp+273}K", expn, ppmRange])

# infos.append([f"Hanecho_{30+273}K_bis", 47, ppmRange])
# infos.append([f"Hanecho_{10+273}K_bis", 88, ppmRange])
# infos.append([f"Hanecho_{0+273}K_bis", 98, ppmRange])
# infos.append([f"Hanecho_{-10+273}K_bis", 109, ppmRange])
# infos.append([f"Hanecho_{0+273}K_bis", 117, ppmRange])
# infos.append([f"Hanecho_{10+273}K_bis", 126, ppmRange])
# infos.append([f"Hanecho_{20+273}K_bis", 136, ppmRange])


# savepath_especifico = "7Li/"
# # info: muestra, expn, ppmRange
# nucleo = "7Li"
# ppmRange = [250, -250]
# infos = [
#     ['zg_293K', 33, ppmRange], # 20 degC
#     ['zgdec_293K', 34, ppmRange],
#     ['zg_293K_longerAcq', 35, ppmRange],

#     ['zg_303K', 43, ppmRange], # 30 degC
#     ['zgdec_303K', 44, ppmRange],
#     ['zg_303K_longerAcq', 45, ppmRange],

#     ['zg_313K', 53, ppmRange], # 40 degC
#     ['zgdec_313K', 54, ppmRange],
#     ['zg_313K_longerAcq', 55, ppmRange],
#     ['zg_313K_evenLongerAcq', 56, ppmRange],

#     ['zg_303K', 63, ppmRange], # 30 degC
#     ['zgdec_303K', 64, ppmRange],
#     ['zg_303K_longerAcq', 65, ppmRange],

#     ['zg_293K', 73, ppmRange], # 20 degC
#     ['zgdec_293K', 74, ppmRange],
#     ['zg_293K_longerAcq', 75, ppmRange],

#     ['zg_283K', 83, ppmRange], # 10 degC
#     ['zgdec_283K', 84, ppmRange],
#     ['zg_283K', 85, ppmRange],

#     ['zg_273K', 93, ppmRange], # 0 degC
#     ['zgdec_273K', 94, ppmRange],
#     ['zg_273K', 95, ppmRange],

#     ['zg_263K', 103, ppmRange], # -10 degC
#     ['zgdec_263K', 104, ppmRange],
#     ['zg_263K', 105, ppmRange],

#     ['zg_273K', 113, ppmRange], # 0 degC
#     ['zgdec_273K', 114, ppmRange],
#     ['zg_273K', 115, ppmRange],

#     ['zg_283K', 123, ppmRange], # 10 degC
#     ['zgdec_283K', 124, ppmRange],
#     ['zg_283K', 125, ppmRange],

#     ['zg_293K', 133, ppmRange], # 20 degC
#     ['zgdec_293K', 134, ppmRange],
#     ['zg_293K', 135, ppmRange],
# ]
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


info_list = []
fig, ax = plt.subplots(num=1, nrows=1, ncols=1)  # create figure & 1 axis
for info in infos:
    muestra, expn, ppmRange= info
    save = True
    muestra = ""
    print("WARNING!!!! Sobreescribiendo muestra a vacio para que no guarde el nombre de la muestra en el archivo de datos")


    # extraigo:
    datos = DatosProcesados(f'{path_local}{path_bruker}{expn}/')
    # datos = DatosProcesados(f'{path_local}')
    if ppmRange is not None:
        datos.espectro.ppmSelect(ppmRange)

    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis

    if muestra.startswith('CP_'):
        contactTime = datos.acqus.dic['P'][15]
        info_list.append([expn, contactTime])


    # grafico para guardar:
    ax.plot(ppmAxis, re, linewidth=2)
    ax.axhline(0, color='k', ls='--')
    ax.set_title(muestra)
    ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
    # ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    if muestra != "":
        muestra = f"_{muestra}"
    if save:
        savepath = f"{savepath_local}{savepath_especifico}"
        # guardo:
        header = "ppmAxis\t real \t imag"
        dataexport = np.array([ppmAxis, re, im]).T
        filename = f'{savepath}/{nucleo}_Nexp{expn}{muestra}.dat'
        np.savetxt(filename, dataexport, header=header)

        # filename = f'{savepath}/{nucleo}_{muestra}.png'
        # fig.savefig(filename)   # save the figure to file
        # # plt.close(fig)    # close the figure window

plt.show()

np.savetxt(f"{savepath_local}{savepath_especifico}contact_times.dat", np.array(info_list), header="expn\t contactTime_us")

# vfit=VoigtFit(ppmAxis, 
#               re, 
#               Npicos=4,
#               ajustar=True,
#               center=[0, 0,-1,-3]
#               )
# fig = vfit.plot_ajuste()
# fig.gca().set_xlim([8,-8])










