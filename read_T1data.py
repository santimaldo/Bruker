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


expn = 6
path = f"S:/NMRdata/2024_Carbones_Fran/2024-04-17_carbones_Fran/{expn}/"
muestra = "1H_0pc_T1"
savepath = "G:/Otros ordenadores/Oficina/Posdoc/CarbonesFran/Datos_Bruker/T1/"

save = True
plotRange = [10, -5]
# rango de integracion
ppmRanges = [[7.5, -2.5],
             [5, 4],
             [2, 1.4],
             [7, 8],
             [9, 10],
             [-2, -4.5]]

datos = DatosProcesadosT1(path)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

# obtengo el espectro para el ultimo variable delay:
re = datos.espectro.real[-1]
im = datos.espectro.imag[-1]

fig_spec, ax_spec = plt.subplots(num=17856)
ax_spec.plot(ppmAxis, re)

colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
Signals = []
ii = -1
# fig, axs = plt.subplots(2, 2)
for ppmRange in ppmRanges:
    ii += 1
    color = colors[ii]
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
    ax_spec.axhline(0, color='k')

    tau, signal = datos.get_T1data(ppmRange)
    tau_fit, signal_fit, residuals = datos.T1fit()
    Signals.append(signal)

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(muestra)
    # -------------
    axs[0, 0].plot(tau, signal, 'o', color=color)
    axs[0, 0].plot(tau_fit, signal_fit, '-', color=color)
    text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
                   multialignment="left")
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    # -------------
    axs[1, 0].plot(tau, residuals, 'o', color=color)
    axs[1, 0].axhline(0, color='k', linestyle='--')
    axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
    # -------------
    axs[0, 1].plot(tau, signal, 'o', color=color)
    axs[0, 1].plot(tau_fit, signal_fit, '-', color=color)
    axs[0, 1].set(xlabel=r'$\tau$ [ms]')
    axs[0, 1].set_xscale('log')
    # -------------
    axs[1, 1].plot(tau, residuals, 'o', color=color)
    axs[1, 1].axhline(0, color='k', linestyle='--')
    axs[1, 1].set(xlabel=r'$\tau$ [ms]')
    axs[1, 1].set_xscale('log')

    for ax in axs.flat:
        ax.label_outer()

# guardo data:
if save:
    filename = f'{savepath}/{muestra}_T1.png'
    fig.savefig(filename)   # save the figure to file

    Signals = np.array(Signals).T
    tau = tau.reshape(tau.size, 1)
    T1data = np.hstack((tau, Signals))
    header = "tau [s]\t"
    for ppmRange in ppmRanges:
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)
