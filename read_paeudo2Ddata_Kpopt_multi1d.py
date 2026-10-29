# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

from paths import NMRDATA, PROJECTS
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

def ppm_correction(ppmAxis, bsms, alpha, beta, ppm_ref):
    """
    Tengo que haber calibrado el eje ppm vs bsms para el mismo SR.
    asi obtengo ppm = alpha * bsms + beta EN LA REFERENCIA
    luego calculo el delta para que la refe siempre caiga en ppm_ref,
    y eso es lo que tengo que correr el ppmAxis
    """
    Delta = (alpha*bsms+beta) - ppm_ref
    return ppmAxis-Delta




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
path  =rf"{NMRDATA}\400dnp\2026-01-22_InSitu/"
# directorio de guradado
savepath= rf"{PROJECTS}\LiMetal\Bruker\analysis\2026_01_InSitu/"

expns = np.arange(24, 45)
bsms_list = np.array([-5750, -9999, -9000, -8750, -7000, -6000, -5000, -4500, 9999, -5500, -8000, -4000, -3000, -2000, 0, 4000, 5000, 1000, 2000, 500, -1000])
savepath += "kpopt/"
muestra = "LiCuFoil_celgard_NMC"
# expns = [1]
# bsms_list = np.array([0])
# muestra = "000_LiCuFoil_celgard_NMC_uwOFF"


save = True
plotRange = [700,-200]
# rango de integracion
ppmRanges = []
window_width = 200  # ancho de cada ventana desde el mínimo local
# parametros de la calibracion hecha previamente
alpha = -1.168e-2
beta = 269.946
ppm_ref = 270 


#=====================================================================
# 2D experiments
#=====================================================================
# fig_popt, ax_popt = plt.subplots(num=382910)
fig_spec, ax_spec = plt.subplots(num=17856)   
fig1, ax1 = plt.subplots()

for jj, expn in enumerate(expns):
   
    bsms = bsms_list[jj]
    datos = DatosProcesados(f'{path}/{expn}/',
                                read_pp = False,
                                nomrmalize_NS=True)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    ppmAxis = ppm_correction(ppmAxis, bsms, alpha, beta, ppm_ref)
    spec1d = datos.espectro.real
    speci1d = datos.espectro.imag

    ax_spec.plot(ppmAxis, spec1d)

    np.savetxt(f"{savepath}/{muestra}_bsms_{int(bsms):04d}.dat", np.array([ppmAxis, spec1d]).T, header="ppmAxis\treal")
    maximo = ppmAxis[spec1d==np.max(spec1d)]
    # maximo = datos.acqus.SFO1



    ax1.scatter(bsms, maximo)
np.savetxt(f"{savepath}/bsms_list", bsms_list)

        
