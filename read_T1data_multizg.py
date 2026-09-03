# -*- coding: utf-8 -*-
"""
Iterar y guardar espectros 1D en un rango de ppm específico
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from Datos import *
# Guardar archivos
save = True
# ------------------ Configuración ------------------
# ########### R13 - NMC-Cu PC protocol
# directorio de datos
expns = np.arange(35, 920, 10)
path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-31_ATMC1_Rui-R13_NMC-Cu_PC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R13_PC/T1data/"
muestra = "Cell_R13_electrolyte"
# rango de integracion
ppmRange = [30, -30]
# Rango de ppm que quieres seleccionar
plotRange = [50, -50]
# # - - - - - - - - - - - - - - - - - - - -
# # ########### R14 - NMC-Cu PC protocol
# # directorio de datos
# expns = np.arange(35, 920, 10)
# path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-31_ATMC1_Rui-R13_NMC-Cu_PC/"
# # directorio de guradado
# savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R13_PC/T1data/"
# muestra = "Cell_R13_electrolyte"
# # rango de integracion
# ppmRange = [30, -30]
# # Rango de ppm que quieres seleccionar
# plotRange = [50, -50]
# # - - - - - - - - - - - - - - - - - - - -

save = True
os.makedirs(savepath, exist_ok=True)
T1s = []
figspec, axspec = plt.subplots(figsize=(8,5))
fig, ax = plt.subplots(figsize=(8,5))
for expn in expns:
    # ------------------ Carga de datos ------------------
    datos = DatosProcesadosT1(f"{path}/{expn}/", vdlist_path="lists/vd/sm2974-T1")
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    last_spec = datos.espectro.real[-1]
    tau, signal = datos.get_T1data(ppmRange=ppmRange)
    tau = tau # paso a segundos

    tau_fit, signal_fit, residuals = datos.T1fit()
    T1s.append(datos.T1params[1])
    
    ax.plot(tau, signal, 'o', label=f"Expn {expn}")
    ax.plot(tau_fit, signal_fit, '-')
    axspec.plot(ppmAxis, last_spec, label=f"Expn {expn}")
    T1datafile = f"{savepath}/T1data_{expn:03d}.dat"
    np.savetxt(T1datafile, np.column_stack((tau, signal)), header="tau [s]\tSignal")
ax.legend()
ax.set_xscale('log')
# ------------------ Graficar todos los espectros ------------------

#%%
fig, ax = plt.subplots(figsize=(8,5))
ax.plot(expns, T1s, 'o-')
for x, y, label in zip(expns, T1s, expns):
    ax.annotate(f"{label}", (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
# %%
