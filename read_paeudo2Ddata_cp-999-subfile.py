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
##----------------------------------------------------------
# # directorio de datos: LiOH (anhydrous)
expn = 44
Npopts = 4
parameter = "P 15"
path  =rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/400dnp/3.2mm-Santi-IMECdendrites-2025-04-28/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\IMEC\DNP\2025-04-28_CP\LiOH_CP-vs-ContactTime/"
sample = "LiOHanhyd"
savefile = f"{sample}_cpHtoLi_uW-OFF_P15"
save = False
plotRange = [-230,-310]
# rango de integracion
ppmRanges = [[-240, -300]]
####################### end Select experiment #######################

#=====================================================================
# 2D experiments
#=====================================================================
popt_subfiles = np.arange(999, 999-Npopts, -1) # 999, 998, ..., 1000-n.
popt_parlist = np.array([])

Signals = np.array([])
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856)
fig_spec_n, ax_spec_n = plt.subplots(num=17857)

datos = DatosProcesados2D(f'{path}/{expn}999/',
                            read_pp = False)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

for popt_subfile in popt_subfiles:
    popt_par = extract_popt_parameters(popt_subfile, parameter, path=path, expn=expn)
    # complete the par list with zeros if the size is not equal to the number of spectra
    popt_parlist = np.append(popt_parlist, popt_par)
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
        ax_spec.axvspan(r1, r2, alpha=0.15, color=color)
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
#%%   
# guardo data:
if save:
    filename = f'{savepath}/{sample}_T1.png'
    fig.savefig(filename)   # save the figure to file

    Signals = np.array(Signals).T
    tau = tau.reshape(tau.size, 1)
    T1data = np.hstack((tau, Signals))
    header = "tau [s]\t"
    for ppmRange in ppmRanges:
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{sample}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{sample}_ultimoEspectro.dat", data)


# fig,ax = plt.subplots(num=1785731)
# ax.plot(popt_parlist[:ppm_of_max.size], ppm_of_max, 'o-')
# ax.axhline(ppm_of_max_in_equilibrium, color='k', linestyle='--')
# ax.set_xlabel('Contact Time [s]')
# ax.set_ylabel(r'something')
# ax.set_xscale('log')
# %%
