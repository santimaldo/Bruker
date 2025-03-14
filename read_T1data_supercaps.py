# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *


# Data paths
expn = 201
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-03-10_insitu-LiTFSIaq-supercap\{expn}/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\sr/"
muestra = f"Nexp{expn}"
save = True
plotRange = [4, -8]

# rango de integracion
ppmRanges = [[3,-6]
             ]

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
figures = []
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
    figures.append(fig)

# guardo data:
if save:
    for ii, fig in enumerate(figures):
        filename = f'{savepath}/{muestra}_T1_range{ii}.png'
        fig.savefig(filename)   # save the figure to file

    Signals = np.array(Signals).T
    tau = tau.reshape(tau.size, 1)
    T1data = np.hstack((tau, Signals))
    header = "tau [s]\t"
    for ii, ppmRange in enumerate(ppmRanges):
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{muestra}_lastSpectrum.dat", data)
    fig_spec.savefig(f"{savepath}/{muestra}_lastSpectrum_withRanges.png")
