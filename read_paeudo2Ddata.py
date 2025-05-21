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
expn_before = 81
expn_pseudo2d = 82
expn_after = 69

path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-02_LiTFSI1M-aq_CA-cycles/"
muestra = "19F_chronoamperometry_0Vto1V"

save = False
plotRange = [5, -10]
# rango de integracion
ppmRanges = [[4, -9],
             [4, -0.5],
             [-0.5, -9]            
            ]

#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
datos = DatosProcesados(f'{path}/{expn_before}/')
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
re = datos.espectro.real
ppmAxis = datos.espectro.ppmAxis
# find the ppm of the maximum in the range < ppm_trershold ppm
ppm_treshold = -1.9
re_in_ROI = re[ppmAxis < ppm_treshold] # re in Region Of Interest
ppmAxis_in_ROI = ppmAxis[ppmAxis < ppm_treshold]
max_index = np.argmax(re_in_ROI)
ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]


#=====================================================================
# Ajuste de espectro antes del experimento 2D
#=====================================================================
# rango de integracion
datos = DatosProcesados2D(f'{path}/{expn_pseudo2d}/')
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

try:
    tau = datos.get_vdlist()/1000  # en segundos
except:
    tau = datos.acqus.D1 * np.arange(1, datos.acqu2s.TD+1)  # Convert to seconds


# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)
ax_spec.axvline(x=ppm_treshold, color='k', linestyle='--')
ppm_of_max = []
for ii in range(datos.espectro.size[0]):
    re = datos.espectro.real[ii]
    im = datos.espectro.imag[ii]
    ax_spec.plot(ppmAxis, re)
    # find the ppm of the maximum in the range < ppm_treshold ppm
    re_in_ROI = re[ppmAxis < ppm_treshold] # re in Region Of Interest
    ppmAxis_in_ROI = ppmAxis[ppmAxis < ppm_treshold]
    max_index = np.argmax(re_in_ROI)
    ppm_of_max_in_equilibrium = ppmAxis_in_ROI[max_index]
    ppm_of_max.append(ppm_of_max_in_equilibrium)
ppm_of_max = np.array(ppm_of_max)


fig,ax = plt.subplots(num=1785731)
ax.plot(tau, ppm_of_max[:tau.size], 'o-')
ax.axhline(ppm_of_max_in_equilibrium, color='k', linestyle='--')
ax.set_xlabel('variable delay [s]')
ax.set_ylabel(r'in-pore $\Delta\delta$ [ppm]')
ax.set_xscale('log')


colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
Signals = []
ii = -1
for ppmRange in ppmRanges:
    ii += 1
    print(ii)
    color = colors[ii]
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
    ax_spec.axhline(0, color='k')

    signal = datos.Integrar(ppmRange=ppmRange)
    #tau_fit, signal_fit, residuals = datos.T1fit()
    Signals.append(signal)

    fig, ax = plt.subplots(num=(ii+1)*382910)
    fig.suptitle(muestra)
    # -------------
    # nd = tau.size
    # ax.plot(tau, signal[:int(nd/2)], 'o', color=color, label="Rising Edge")
    #   ax.plot(tau, signal[int(nd/2):], 'o', color=color, alpha=0.3, label="Falling Edge")
    # ax.plot(tau, signal, 'o', color=color, label="Rising Edge")
    ax.plot(tau, signal[:tau.size], 'o', color=color, label="Rising Edge")
    ax.set_xscale('log')
    #axs[0, 0].plot(tau_fit, signal_fit, '-', color=color)
    #text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
    #axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
    #               multialignment="left")
    #axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    
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
