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
import re

################## Functions ###########################
def find_popt_subfolders(expn, path):
    """
    Find all subfolders that start with the given 
    experiment number followed by 99n:
    999, 998, ..., 1000-n.
    """
    experiment_numbers = []
    for folder in os.listdir(path):
        if folder.startswith(f"{expn}"):
            try:
                suffix = int(folder[len(f"{expn}"):])
                experiment_numbers.append(suffix)
            except ValueError:
                continue
    return np.array(sorted(experiment_numbers, reverse=True))

import numpy as np

def extract_popt_parameters(parameter="p15", group=0, path=".", expn=1):
    """
    Extracts and concatenates parameter arrays from a popt.protocol file.

    Each matching line (Step:) filtered by group and parameter is used to
    generate a np.linspace array from 'desde' to 'hasta' with 'num_puntos' points.
    All such arrays are concatenated into one 1D numpy array.

    Parameters:
        subfile (str): suffix of the file to read (path/expn/popt.protocol.subfile)
        parameter (str): parameter name to filter lines (default 'p15')
        group (int): group number to filter lines (default 0)
        path (str): base directory (default '.')
        expn (int): experiment subfolder (default 1)

    Returns:
        np.ndarray: concatenated 1D array of parameter values.
    """
    filename = f"{path}/{expn}/popt.array"
    arrays = []

    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith("Step:"):
                continue
            parts = line.strip().split()
            try:
                line_group = int(parts[1])
                line_param = parts[2]
                desde = float(parts[4])
                hasta = float(parts[5])
                num_puntos = int(parts[6])
            except (ValueError, IndexError):
                print(f"Skipping line due to parsing error: {line.strip()}")
                continue

            if line_group == group and line_param == parameter:
                arr = np.linspace(desde, hasta, num_puntos)
                arrays.append(arr)

    if arrays:
        return np.concatenate(arrays)
    else:
        return np.array([])

#%%                
    return np.array(data)
################## end Functions #######################


################## Select experiment #######################
# # directorio de datos
# expn = 64
# Npopts = 3
# parameter = "D 1"
# path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\IMECeLi_CP-vs-D1_P15-01ms/"
# sample = "IMECdendrites_KBr"
# savefile = f"{sample}_cpLitoH_uW-ON_D1"
# save = False
# plotRange = [240, 200]
# # rango de integracion
# ppmRanges = [[233, 210]]
##----------------------------------------------------------
# directorio de datos
# expn = 65
# Npopts = 3
# parameter = "D 1"
# path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\IMECeLi_CP-vs-D1_P15-10ms/"
# sample = "IMECdendrites_KBr"
# savefile = f"{sample}_cpLitoH_uW-ON_D1"
# save = False
# plotRange = [240, 200]
# # rango de integracion
# ppmRanges = [[233, 210]]
# ##----------------------------------------------------------
# # # directorio de datos: LiOH (anhydrous)
# expn = 44
# Npopts = 4
# parameter = "P 15"
# path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\LiOH_CP-vs-ContactTime/"
# sample = "LiOHanhyd"
# savefile = f"{sample}_cpHtoLi_uW-OFF_P15"
# save = False
# plotRange = [-230,-310]
# # rango de integracion
# ppmRanges = [[-240, -300]]
##----------------------------------------------------------
# # directorio de datos: LiOH (anhydrous)
expn = 21
Npopts = 1
parameter = "p15"
path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-06-17_3.2mm_IMECdendrites/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\analysis\2025_06_IMEC\CP_7Li-1H_uW-ON/"
sample = "IMECeLi"
savefile = f"{sample}_cpHtoLi_uW-ON_P15"
save = False
plotRange = [50,-50]
# rango de integracion
ppmRanges = [[10, -10]]
####################### end Select experiment #######################

#=====================================================================
# 2D experiments
#=====================================================================

Signals = np.array([])
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)
fig_spec_n, ax_spec_n = plt.subplots(num=17857)

datos = DatosProcesados2D(f'{path}/{expn}999/',
                            read_pp = False)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real


popt_parlist = extract_popt_parameters(parameter=parameter, group=0, path=path, expn=expn)
# complete the par list with zeros if the size is not equal to the number of spectra
popt_parlist = np.append(popt_parlist, np.zeros(spec.shape[0] - popt_parlist.size))
# Plot the 1D spectra
for kk in range(spec.shape[0]):
    ax_spec.plot(ppmAxis, spec[kk,:])            
    if popt_parlist[kk] != 0:
        data = np.array([ppmAxis, spec[kk,:]]).T
        np.savetxt(f"{savepath}/{savefile}_{popt_parlist[kk]:0>4.1f}_s.dat",
                    data, header="ppmAxis\t real")

    ###### start integrating the spectra
    colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
    ii = -1
    for ppmRange in ppmRanges:
        ii += 1
        color = colors[ii]
        ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
        r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
        ax_spec.axvline(r1, color='k', linestyle='--')
        ax_spec.axvline(r2, color='k', linestyle='--')
        ax_spec.axhline(0, color='k')

        signal = datos.Integrar(ppmRange=ppmRange) 
        #tau_fit, signal_fit, residuals = datos.T1fit()
        Signals = signal
    


popt_parlist = popt_parlist[Signals!=0] # remove zeros generated by stopping the experiment
Signals = Signals[Signals != 0] # remove zeros generated by stopping the experiment


fig_popt, ax_popt = plt.subplots(num=382910)
ax_popt.plot(popt_parlist/1000, Signals, 'o')#, color=color, label="Rising Edge")
ax_popt.xaxis.set_major_locator(plt.MultipleLocator(1))
ax_popt.set_xlabel("Contact Time [ms]")
ax_popt.set_ylabel("Integral [a.u.]")
ax_popt.set_title(r"$^1$H $\rightarrow$ $^7$Li CP - LiOH")
