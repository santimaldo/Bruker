# -*- coding: utf-8 -*-
"""
Adapted on Sat Nov 8 2025

Reemplaza la integración por ajuste de dos picos pseudo-Voigt.
Guarda los resultados en un DataFrame ordenado por contact time.
@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Datos import *
from VoigtFit import PseudoVoigtFit
import re

#=====================================================================
# Directorios y parámetros
#=====================================================================

expns = np.concatenate([np.arange(60, 65), np.arange(66, 71)])
expns = [69]
path  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/"
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-02_LiTFSI1M-aq_CA-cycles/"
muestra = ""
save = False

plotRange = [100, -100]
ppmRange  = [20, -15]   # rango de ajuste
show_individual_fits = True  # para graficar cada ajuste individual

#=====================================================================
# Inicialización
#=====================================================================

results = []   # lista donde acumularé resultados fila por fila

#=====================================================================
# Bucle sobre experimentos
#=====================================================================

fig_spec, ax_spec = plt.subplots(num=17856)

for jj, expn in enumerate(expns):

    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados(f'{path}/{expn}/', read_pp=False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    contact_time = datos.acqus.dic["P"][15]

    # Graficar espectro completo
    ax_spec.plot(ppmAxis, spec, lw=1)
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax_spec.axvspan(np.min(ppmRange), np.max(ppmRange), alpha=0.15, color='gray')
    ax_spec.axhline(0, color='k')

    # Seleccionar solo la región de ajuste
    mask = (ppmAxis <= np.max(ppmRange)) & (ppmAxis >= np.min(ppmRange))
    xdata = ppmAxis[mask]
    ydata = spec[mask]

    # guess from topspin:
    # amplitude = [011102470.3/(datos.acqus.NS*datos.acqus.RG),
    #              444291004.2/(datos.acqus.NS*datos.acqus.RG)]
    # center = [5.251, 0.00] # ppm
    # fwhm = [3013.5789 / datos.acqus.SFO1,
    #         3230.6543 / datos.acqus.SFO1]
    # fraction = [0.3296, 0.4225]  # Lorentzian fraction
    amplitude = [466021696.3/(datos.acqus.NS*datos.acqus.RG),
                 604812667/(datos.acqus.NS*datos.acqus.RG)]
    center = [4, -0.65] # ppm
    fwhm = [2602.5823 / datos.acqus.SFO1,
            3220.4293 / datos.acqus.SFO1]
    fraction = [0.6742, 0.0772]  # Lorentzian fraction
    vfit = PseudoVoigtFit(xdata,
                    ydata,
                    Npicos=2,
                    ajustar=True,
                    amplitude=amplitude,
                    center=center,
                    fwhm=fwhm,
                    fraction=fraction,
                    fijar=["center", "fwhm", "fraction"]
                    )

    # Extraer amplitudes y errores
    m1_area = vfit.params["m1_amplitude"].value
    m2_area = vfit.params["m2_amplitude"].value
    m1_area_err  = vfit.params["m1_amplitude"].stderr
    m2_area_err  = vfit.params["m2_amplitude"].stderr
    m1_center = vfit.params["m1_center"].value
    m2_center = vfit.params["m2_center"].value
    m1_center_err = vfit.params["m1_center"].stderr
    m2_center_err = vfit.params["m2_center"].stderr


    # Guardar resultados en la lista
    results.append({
        "expn": expn,
        "contact_time": contact_time,
        "m1_area": m1_area,
        "m1_area_err": m1_area_err,
        "m2_area": m2_area,
        "m2_area_err": m2_area_err,
        "m1_center": m1_center,
        "m1_center_err": m1_center_err,
        "m2_center": m2_center,
        "m2_center_err": m2_center_err
    })

    # Mostrar ajuste individual
    if show_individual_fits:
        fig_fit = vfit.plot_ajuste(xlabel=r"$^1$H Chemical Shif [ppm]",
                                   ylabel="Intensity (a.u.)",
                                   reverse_xaxis=True)
        fig_fit.gca().set_title(f"Exp {expn} — Contact time = {contact_time:.1f} μs")


#=====================================================================
# Crear DataFrame y ordenarlo
#=====================================================================

df_results = pd.DataFrame(results)
df_results.sort_values(by="contact_time", inplace=True, ignore_index=True)

print(df_results)

#%%=====================================================================
# Gráfico de resultados
#=====================================================================
colors = ['r', 'b']
fig_area, ax_area = plt.subplots(num=382910)
ax_area.errorbar(df_results["contact_time"], df_results["m1_area"], yerr=1.96*df_results["m1_area_err"], fmt='o-', label=f'Peak: {center[0]}', color=colors[0])
ax_area.errorbar(df_results["contact_time"], df_results["m2_area"], yerr=1.96*df_results["m2_area_err"], fmt='o-', label=f'Peak: {center[1]}', color=colors[1])
ax_area.set_xlabel("Contact time")
ax_area.set_ylabel("Amplitude (a.u.)")
ax_area.legend()
ax_area.grid(True)



#=====================================================================
# Guardado de resultados
#=====================================================================

if save:
    fig_area.savefig(f"{savepath}/{muestra}_Areas_vs_contactTime.png")
    df_results.to_csv(f"{savepath}/{muestra}_VoigtFit_results.csv", index=False)

# %%
