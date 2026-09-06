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
# # ########### R13 - NMC-Cu PC protocol
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
# # # - - - - - - - - - - - - - - - - - - - -
# ########### R14 - NMC-Cu CC protocol
# directorio de datos
# expns = np.arange(35, 800, 10)
expns = np.concatenate( (np.arange(35, 890, 10), np.arange(1005, 1442, 10)) )
path = rf"c:\Users\Santi\Documents\NMRdata\300neo\2026-09-02_ATMC1_Rui-R14_NMC-Cu_CC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-09_R14_CC/T1data/"
muestra = "Cell_R14_electrolyte"
# rango de integracion
ppmRange = [30, -30]
# Rango de ppm que quieres seleccionar
plotRange = [50, -50]
# - - - - - - - - - - - - - - - - - - - -

save = True
os.makedirs(savepath, exist_ok=True)
os.makedirs(savepath+"figs/", exist_ok=True)
T1s = []
times = []
figspec, axspec = plt.subplots(figsize=(8,5))
fig, ax = plt.subplots(figsize=(8,5))
for i, expn in enumerate(expns):
    # ------------------ Carga de datos ------------------
    datos = DatosProcesadosT1(f"{path}/{expn}/", vdlist_path="lists/vd/sm2974-T1")
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    last_spec = datos.espectro.real[-1]
    tau, signal = datos.get_T1data(ppmRange=ppmRange)
    tau = tau # paso a segundos
    times.append(datos.acqus.dic["DATE"])
    if i==0:
        t0 = datos.acqus.dic["DATE"]
    time_rel_h = (datos.acqus.dic["DATE"]-t0)/3600

    tau_fit, signal_fit, residuals = datos.T1fit()
    T1s.append(datos.T1params[1])
    T1datafile = f"{savepath}/T1data_{expn:03d}.dat"
    #np.savetxt(T1datafile, np.column_stack((tau, signal)), header="tau [s]\tSignal")

    
    color = plt.cm.viridis(i / max(len(expns) - 1, 1))
    axspec.plot(ppmAxis, last_spec, color=color, label=f"Expn {expn}")

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(f"{muestra} - time {time_rel_h:.1f} h")
    # -------------
    axs[0, 0].plot(tau, signal, 'ko')
    axs[0, 0].plot(tau_fit, signal_fit, 'r-')
    text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
                   multialignment="left")
    axs[0, 0].set(xlabel=r'time [ms]', ylabel=r'$S$')
    # -------------
    axs[1, 0].plot(tau, residuals, 'ko')
    axs[1, 0].axhline(0, color='r', linestyle='--')
    axs[1, 0].set(xlabel=r'time [ms]', ylabel=r'Residuals')
    # -------------
    axs[0, 1].plot(tau, signal, 'ko')
    axs[0, 1].plot(tau_fit, signal_fit, 'r-')
    axs[0, 1].set(xlabel=r'time [ms]')
    axs[0, 1].set_xscale('log')
    # -------------
    axs[1, 1].plot(tau, residuals, 'ko')
    axs[1, 1].axhline(0, color='r', linestyle='--')
    axs[1, 1].set(xlabel=r'time [ms]')
    axs[1, 1].set_xscale('log')
    fig.savefig(f"{savepath}/figs/T1fit_{expn:03d}.png", dpi=300)
    plt.close(fig)


# ------------------ Graficar todos los espectros ------------------
savedata = np.array([times, T1s]).T
np.savetxt(f"{savepath}/T1s_vs_time.txt", savedata, header="Time\tT1 [ms]")
#%%
fig, ax = plt.subplots(figsize=(8,5))
ax.plot(expns, T1s, 'o-')
for x, y, label in zip(expns, T1s, expns):
    ax.annotate(f"{label}", (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
# %%
