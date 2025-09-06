# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Autophase import autophase
import scipy.integrate as integrate
import re
import shutil
import os

################## Functions ###########################

def read_bsms_field(path_archivo):
    """
    Lee un archivo con formato específico y devuelve el array 'x' como lista de floats.

    Parámetros:
    path_archivo (str): Ruta al archivo de entrada.

    Retorna:
    list[float]: Lista con los valores del array x.
    """
    with open(path_archivo+'Klog_opt', 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                x = eval(contenido)
                return x
    raise ValueError("No se encontro una linea que empiece con 'x='.")


import numpy as np
from scipy import integrate

def integrate_around_peak(ppmAxis, spectrum, window_width, range_of_max=None):
    """
    Integrates a spectrum around the maximum peak within a specified window width.

    Parameters:
    - ppmAxis (np.array): ppm axis (can be decreasing or increasing)
    - spectrum (np.array): corresponding real spectrum
    - window_width (float): integration window in ppm
    - range_of_max (tuple, optional): (ppm_min, ppm_max) to search for the maximum. 
                                      If None, use the full spectrum.

    Returns:
    - float: integral value in the range centered at the maximum
    """
    # Mask to restrict search for the maximum if requested
    if range_of_max is not None:
        mask_max = (ppmAxis >= min(range_of_max)) & (ppmAxis <= max(range_of_max))
        if not np.any(mask_max):
            raise ValueError("No points found in the specified range_of_max.")
        idx_max = np.argmax(np.abs(spectrum[mask_max]))
        max_index = np.where(mask_max)[0][idx_max]
    else:
        max_index = np.argmax(np.abs(spectrum))

    center_ppm = ppmAxis[max_index]

    # Define integration window
    r1 = center_ppm - window_width / 2
    r2 = center_ppm + window_width / 2

    # Handle both decreasing and increasing ppm axis
    lower, upper = sorted([r1, r2])
    mask = (ppmAxis >= lower) & (ppmAxis <= upper)

    # Bruker convention: positive peaks down, hence negative integral
    return -integrate.simpson(spectrum[mask], x=ppmAxis[mask])



#%%                
    return np.array(data)
################## end Functions #######################



# directorio de datos
expn = 16
absolute = False
autoph = True 
path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-07-31_InSitu/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker/analysis/"
muestra = ""
save = False
plotRange = [1000,-400]
# rango de integracion
ppmRanges = [[800, 150],
             [150,-200],
            #[-0.5, -9]            
            ]
window_width = 50  # ancho de cada ventana desde el mínimo local


#=====================================================================
# 2D experiments
#=====================================================================
fig_spec, ax_spec = plt.subplots(num=17856)
fig_popt, ax_popt = plt.subplots(num=382910)





path_2D = f"{path}/{expn}/"
acqus_file = os.path.join(path_2D, "acqus")

bsms_field = np.array(read_bsms_field(path_2D))
Signals = np.zeros([bsms_field.size, len(ppmRanges)])

for kk in range(bsms_field.size):
    # copy acqus file. Esto es necesario para que DatosProcesados lea bien los datos
    carpeta_destino = os.path.join(path, f"{expn}{kk+1:02d}/")
    shutil.copy(acqus_file, carpeta_destino)

    datos = DatosProcesados(f'{path}/{expn}{kk+1:02d}/',
                                read_pp = False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec1d = datos.espectro.real
    speci1d = datos.espectro.imag

    # if autoph:
    #     spec1d = ng.proc_autophase.autops(spec1d+1j*speci1d, "acme") 
    #     spec[kk,:]  = ng.process.proc_bl.cbf(spec1d.real, last=100)
    ax_spec.plot(ppmAxis, spec1d)

    if kk==0:
        spec = np.zeros([bsms_field.size, spec1d.size])
    spec[kk,:] = spec1d

###### start integrating the spectra
colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
for ii, ppmRange in enumerate(ppmRanges):
    color = colors[ii]
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
    ax_spec.axhline(0, color='k')
    # signal = datos.Integrar(ppmRange=ppmRange)
    for kk in range(bsms_field.size):
        Signals[kk, ii] = integrate_around_peak(ppmAxis, spec[kk, :], window_width, range_of_max=ppmRange)
    


for ii, ppmRange in enumerate(ppmRanges):
    ax_popt.plot(bsms_field, Signals[:,ii], 'o')#, color=color, label="Rising Edge")

# if save:
#     for kk in range(Signals.size):
#         np.savetxt(f"{savepath}/{muestra}_bsms_{bsms_field[kk]}.dat",
#                 np.array([ppmAxis, spec[kk,:]]).T,
#                 header="ppmAxis\treal")
