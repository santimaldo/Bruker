# -*- coding: utf-8 -*-
"""
Created on Wen 19/02/2025

@author: Santi

Experimentos pseudo 2D:
chronoamperometrias y espectros en funcion del tiempo 
de espera guardado en vdlist.

Primero ajustamos el espectro antes del experimento 2D.
Luego, los parametros obtenidos en cada espectro son
usados para el guess del espectro siguiente

"""

import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from VoigtFit import VoigtFit
from scipy.signal import find_peaks


# directorio de datos
expn_before = 67
expn_pseudo2d = 68
expn_after = 69

path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/300old/2025-02-07_insitu-sync-start/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-02_LiTFSI1M-aq_CA-cycles/"
muestra = "19F_chronoamperometry_0Vto1V"

save = False
ppmRange = [10, -10]

Npeaks = 3
center_guess = [0, -0.16, -2.7]
sigma_guess = [0.16, 0.5, -0.76]

#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
datos = DatosProcesados2D(f'{path}/{expn_pseudo2d}/')
datos.espectro.ppmSelect(ppmRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real
re = spec[0,:]
ppmAxis = datos.espectro.ppmAxis

vfit=VoigtFit(ppmAxis, 
              re, 
              Npicos=Npeaks,
              ajustar=True,
              center=center_guess,
              sigma=sigma_guess 
              )
fig = vfit.plot_ajuste()
fig.gca().set_xlim(ppmRange)

#=====================================================================
# Ajuste de espectro antes del experimento 2D
#=====================================================================
# rango de integracion
datos = DatosProcesados2D(f'{path}/{expn_pseudo2d}/')
datos.espectro.ppmSelect(ppmRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real
vfit2d = vfit

# grafico todos los espectros juntos
centers = []
centers_error = []
areas = []
for ii in range(datos.espectro.size[0]):
# for ii in range(5):
    #fig_spec, ax_spec = plt.subplots(num=ii)
    re = datos.espectro.real[ii]
    im = datos.espectro.imag[ii]
    #ax_spec.plot(ppmAxis, re)


    center_guess = [vfit2d.params[f'm{jj+1}_center'].value for jj in range(Npeaks)]
    sigma_guess = [vfit2d.params[f'm{jj+1}_sigma'].value for jj in range(Npeaks)]
    vfit2d=VoigtFit(ppmAxis, 
                    re, 
                    Npicos=Npeaks,
                    ajustar=True,
                    fijar=["m1_center", "m2_center", "m1_sigma", "m2_sigma"],
                    center=center_guess,
                    sigma=sigma_guess
                    )
    fig = vfit2d.plot_ajuste()
    fig.gca().set_xlim(ppmRange)
    centers.append([vfit2d.params[f'm{jj+1}_center'].value for jj in range(Npeaks)])
    centers_error.append([vfit2d.params[f'm{jj+1}_center'].stderr for jj in range(Npeaks)])
    areas.append([vfit2d.params[f'm{jj+1}_amplitude'].value for jj in range(Npeaks)])

centers = np.array(centers)
centers_error = np.array(centers_error)
areas = np.array(areas)

# reordeno filtrando los que tienen mucho error.
# creating a mask
# elements_to_delete = abs(centers/centers_error)>1
# centers[elements_to_delete] = np.nan


plt.figure(1111111)
plt.title("Centers")
plt.plot(centers[:,2],'o-')
for ii in range(Npeaks):

    plt.figure(1111112)
    plt.title("Areas")
    plt.plot(areas[:,ii],'o-')
stop

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
    tau = datos.get_vdlist()/1000
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
