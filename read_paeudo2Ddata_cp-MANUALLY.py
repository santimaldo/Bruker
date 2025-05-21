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

def extract_popt_parameters(subfile, parameter, path=".", expn=1):
    """
    Extracts the specified parameter from the popt.protocol file.
    The parameter is specified as a string.

    typical options are:
    "P 15" for CP contact time
    """
    filename = f"{path}/{expn}/popt.protocol.{subfile}"
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        start_reading = False
        column_index = None
        for line in lines:
            if line.startswith("Experiment") and parameter in line:
                headers = re.split(r'\s{2,}', line)
                column_index = headers.index(parameter)
                start_reading = True
                continue
            if start_reading:
                columns = line.split()
                if columns == [] or columns[0].isdigit() == False:
                    break
                
                data.append(float(columns[column_index]))            
    return np.array(data)

def extract_popt_array(expn=None, path=None):
    values = []
    filename = f"{path}/{expn}/popt.array"
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith("Step:"):
                parts = line.strip().split()
                start = float(parts[4])
                stop = float(parts[5])
                points = int(parts[6])
                step_values = np.linspace(start, stop, points)
                values.append(step_values)

    return np.concatenate(values)
################## end Functions #######################



#################### Select Experiment #######################
# parameter = "D 1"
# path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\IMECeLi_CP-vs-D1_P15-10ms/"
# sample = "IMECdendrites_KBr"
# savefile = f"{sample}_cpLitoH_uW-ON_D1"
# parameter_units = "s" # for the savefile
# save = False
# plotRange = [400, 50]

# expns = [65999, 
#          65996]

# popt_parlists =  []
# # first popt_parlist
# popt_parlist = np.zeros(0)
# for popt_subfile in [999, 998, 997]:
#     popt_par = extract_popt_parameters(popt_subfile, parameter, path=path, expn=65)
#     # complete the par list with zeros if the size is not equal to the number of spectra
#     popt_parlist = np.append(popt_parlist, popt_par)
# popt_parlists.append(popt_parlist)
# # second popt_parlist
# popt_parlist =extract_popt_parameters(996, parameter, path=path, expn=65)
# popt_parlists.append(popt_parlist)
##----------------------------------------------------------
parameter = "D 1"
path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\IMECeLi_CP-vs-D1_P15-01ms/"
sample = "IMECdendrites_KBr"
savefile = f"{sample}_cpLitoH_uW-ON_D1"
parameter_units = "s" # for the savefile
save = False
plotRange = [400, 50]

expns = [64999]

popt_parlists =  []
# first popt_parlist
popt_array = extract_popt_array(path=path, expn=64)
popt_parlists.append(popt_array)

#################### end Select Experiment #######################

#=====================================================================
# 2D experiments
#=====================================================================
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)

cmaps = [plt.cm.Blues, plt.cm.Greens, plt.cm.Reds, plt.cm.Purples, plt.cm.Oranges]
for idx, (expn, popt_parlist) in enumerate(zip(expns, popt_parlists)):
    datos = DatosProcesados2D(f'{path}/{expn}/',
                                read_pp = False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    #fill the popt_parlist with zeros if the size is not equal to the number of spectra
    if popt_parlist.size < spec.shape[0]:
        popt_parlist = np.append(popt_parlist, np.zeros(spec.shape[0] - popt_parlist.size))
    elif popt_parlist.size > spec.shape[0]:
        # remove the last elements of the list
        print("WARNING! popt_parlist is larger than the number of spectra")
        popt_parlist = popt_parlist[:spec.shape[0]]
    # Generate a colormap for the current expn
    cmap = cmaps[idx]
    colors = cmap(np.linspace(0.3, 1, spec.shape[0]))  # Avoid extremes (white/black)

    # Plot the 1D spectra with the colormap
    for kk in range(spec.shape[0]):
        ax_spec.plot(ppmAxis, spec[kk, :], color=colors[kk])
        if popt_parlist[kk] != 0:
            data = np.array([ppmAxis, spec[kk,:]]).T
            np.savetxt(f"{savepath}/{savefile}_{popt_parlist[kk]:0>4.1f}_{parameter_units}.dat",
                        data, header="ppmAxis\t real")

ax_spec.set_xlim(plotRange)
ax_spec.set_xlabel("ppm")