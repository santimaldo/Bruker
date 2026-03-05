# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:00:00 2026

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

# =======================
# Parámetros del script
# =======================
expns = np.arange(20, 36)
path  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2026-01-21_3.2mm_LiH/"
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\Debashis\Analysis\2026-01_LiH\CPspec/"
min_contact_time = 0
muestra = ""
save = True
plotRange = [100, -100]

npts = 6   # puntos de la FID a considerar


# =======================
# Inicialización de arrays
# =======================
contact_times = np.zeros(len(expns))
Signals = np.zeros(len(expns))
Errors = np.zeros(len(expns))


# =======================
# Graficar todos los espectros juntos
# =======================
fig_spec, ax_spec = plt.subplots(num=17856)

for jj, expn in enumerate(expns):
    path_exp = f"{path}/{expn}/"
    datos = DatosProcesados(path_exp, read_pp=False)
    datos.set_fid()
    # Selecciono rango de ppm para mostrar espectro
    contact_time = datos.acqus.dic["P"][15]
    
    if contact_time < min_contact_time:
        continue
    
    contact_times[jj] = contact_time
    
    
    # =======================
    # Tomar valor absoluto de los primeros puntos de la FID
    # =======================
    fid0 = np.abs(datos.fid.real + 1j*datos.fid.imag)
    fid0 = fid0[:32]
    signal = fid0[:npts].mean()       # promedio como señal representativa
    std_err = fid0[:npts].std()/np.sqrt(npts)
    Signals[jj] = signal
    Errors[jj] = std_err

    tau = np.arange(fid0.size) * datos.acqus.DW
    ax_spec.plot(tau, fid0, '.-' )
# =======================
# Graficar CP curve
# =======================
fig_popt, ax_popt = plt.subplots(num=382910)
ax_popt.errorbar(contact_times, Signals, yerr=Errors, fmt='o')
ax_popt.set_xlabel("Contact time [s]")
ax_popt.set_ylabel("Signal amplitude [a.u.]")
ax_popt.grid(True)
ax_popt.set_xscale('log')

# =======================
# Guardar datos si corresponde
# =======================
if save:
    CPdata = np.array([contact_times, Signals]).T
    CPdata = CPdata[CPdata[:,1] != 0]
    CPdata = CPdata[CPdata[:,0].argsort()]
    
    header = "contact_time [s]\tSignal [a.u.]\n"
    np.savetxt(f"{savepath}/{muestra}_CP-curve_FID.dat", CPdata, header=header)
