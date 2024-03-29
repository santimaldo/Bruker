# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate


# directorio de datos
# expn = 102
# path = f"S:/Carbones/300MHz/2022-09-01_CarbonesCNEA_M4_LiTFSI/{expn}/"
# # directorio de guradado
# savepath = "S:/Carbones/analisis/2022-09_M4_LiTFSI-G2_100mM/filesT1/"
# muestra = "tmp"  # '1H_M4-NaOH' # "LiTFSI-G2_100mM_bulk", "M4-NaOH", "M4-HF"


# directorio de datos  ### M4 HF muestra SIN MOLER
# expn = 8
# path = f"S:/Carbones/300MHz/2018-09-13_Carbones_Glyma/{expn}/"
# savepath = "S:/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/T1/"
# muestra = "1H_M4-HF_SinMoler"

# # directorio de datos ### M4 HF muestra MOLIDA
# expn = 11
# path = f"S:/CNEA/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/{expn}/"
# savepath = "S:/CNEA/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/T1/"
# muestra = "1H_M4-HF_Molida"

# ------------------------------- Carbnnes CNEA 2023
# directorio de datos ### M4 HF LITFSI 10/10/2023
# expn = 5
# path = f"S:/NMRdata/2018_Carbones_CNEA/2023-10-10_CarbonesCNEA_M4_Diff/{expn}/"
# savepath = "S:/tmp/"
# muestra = "19F_M4-HF-LITFSI_T1"

# expn = 26 #2
# path = f"S:/NMRdata/2018_Carbones_CNEA/2023-10-10_CarbonesCNEA_M4_Diff/{expn}/"
# muestra = "7Li_M4-HF-LITFSI_T1"

# directorio de datos ### M4 HF LITFSI 10/10/2023
# expn = 6
# path = f"S:/NMRdata/2018_Carbones_CNEA/2023-10-11_CarbonesCNEA_M4-HF-LiTF_Diff/{expn}/"
# muestra = "19F_M4-HF-LITF_T1"

# expn = 9
# path = f"S:/NMRdata/2018_Carbones_CNEA/2023-10-11_CarbonesCNEA_M4-HF-LiTF_Diff/{expn}/"
# muestra = "7Li_M4-HF-LITF_T1"


# expn = 1997
# path = f"S:/Doctorado/Carbones/300MHz/2018-09-25_Carbones_Glyma/{expn}/"
# muestra = "7Li_M4-NaOH-LITf_T1"

expn = 4
path = f"S:/NMRdata/2022_Polisulfuros/2023-04-05_Diff_Polisulfuros/{expn}/"
muestra = "Li2S6-TEGDME-DME"


# savepath = "S:/Posdoc/CNEA/Carbones/analisis/2023-10_Difusion/datos/T1/"
savepath = "S:/tmp/"
# ------------------------------- Carbnnes CNEA 2023


# Polisulfuros
# expn = 2
# muestra = "Li2S6-DME-Dia21"
# path = f"S:/NMRdata/2022_Polisulfuros/2023-03-28_Diff_Polisulfuros/{expn}/"
# savepath = "S:/tmp/"
# savepath = "S:/Posdoc/Li-S/Analisis/2023-03_Li2S6-Diff_DME-TEGDME/T1/"

save = True
# rango de integracion
ppmRange = [-1, 1]
# ppmRange = [0, -8]

datos = DatosProcesadosT1(path)
ppmAxis = datos.espectro.ppmAxis


spec = datos.espectro.real

re = datos.espectro.real[-1]
im = datos.espectro.imag[-1]


plt.figure(7532)
plt.plot(ppmAxis, re)
# plt.plot(ppmAxis, im)
plt.xlim(np.max(ppmAxis), np.min(ppmAxis))
r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
plt.axvspan(r1, r2, alpha=0.2, color='red')
plt.axhline(0, color='k')

# %%

tau, signal = datos.get_T1data(ppmRange)
tau_fit, signal_fit, residuals = datos.T1fit()


fig, axs = plt.subplots(2, 2)
fig.suptitle(muestra)
# -------------
axs[0, 0].plot(tau, signal, 'ko')
axs[0, 0].plot(tau_fit, signal_fit, 'r-')
text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
               multialignment="left")
axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
# -------------
axs[1, 0].plot(tau, residuals, 'ko')
axs[1, 0].axhline(0, color='k', linestyle='--')
axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
# -------------
axs[0, 1].plot(tau, signal, 'ko')
axs[0, 1].plot(tau_fit, signal_fit, 'r-')
axs[0, 1].set(xlabel=r'$\tau$ [ms]')
axs[0, 1].set_xscale('log')
# -------------
axs[1, 1].plot(tau, residuals, 'ko')
axs[1, 1].axhline(0, color='k', linestyle='--')
axs[1, 1].set(xlabel=r'$\tau$ [ms]')
axs[1, 1].set_xscale('log')

for ax in axs.flat:
    ax.label_outer()

# guardo data:
if save:
    filename = f'{savepath}/{muestra}_T1.png'
    fig.savefig(filename)   # save the figure to file

    T1data = np.array([tau, signal]).T
    np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)
