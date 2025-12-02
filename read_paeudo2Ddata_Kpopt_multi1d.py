# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
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
################## end Functions #######################





# directorio de datos
absolute = False
autoph = True 
path  =r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-28_InSitu/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Bruker\analysis\2025-11_InSitu\Kpopt\LiCuFoil-celgard-LFP/"

expns = [71, 73, 75, 77, 79]
muestra = "LiCuFoil_celgard_LFP"
save = True
plotRange = [1000,-400]
# rango de integracion
ppmRanges = []
window_width = 200  # ancho de cada ventana desde el mínimo local



#=====================================================================
# 2D experiments
#=====================================================================
fig_popt, ax_popt = plt.subplots(num=382910)

for jj, expn in enumerate(expns):
    fig_spec, ax_spec = plt.subplots(num=17856+jj)    
    path_2D = f"{path}/{expn}/"
    bsms_field = np.array(read_bsms_field(path_2D))

    for kk, bsms in enumerate(bsms_field):
        acqus_file = os.path.join(path_2D, "acqus")
        nexp1d = f"{expn}999{kk+1:02d}"

        # copy acqus file. Esto es necesario para que DatosProcesados lea bien los datos
        carpeta_destino = os.path.join(path, f"{nexp1d}/")
        shutil.copy(acqus_file, carpeta_destino)

        datos = DatosProcesados(f'{path}/{nexp1d}/',
                                   read_pp = False,
                                   nomrmalize_NS=False)
        datos.espectro.ppmSelect(plotRange)
        ppmAxis = datos.espectro.ppmAxis
        spec1d = datos.espectro.real
        speci1d = datos.espectro.imag

        ax_spec.plot(ppmAxis, spec1d)

        np.savetxt(f"{savepath}/{muestra}_bsms_{int(bsms):04d}.dat", np.array([ppmAxis, spec1d]).T, header="ppmAxis\treal")

        
