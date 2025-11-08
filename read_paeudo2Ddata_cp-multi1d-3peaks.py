# -*- coding: utf-8 -*-
"""
Adapted on Sat Nov 8 2025

Reemplaza la integración por ajuste de tres picos pseudo-Voigt.
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
# expns = [69]  # opcional para probar un solo experimento
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
Npicos = 3     # número de picos a ajustar

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

    # Guess inicial para los tres picos (ajustar según tu espectro)
    amplitude = [359451094.1/(datos.acqus.NS*datos.acqus.RG),
                 162708124.8/(datos.acqus.NS*datos.acqus.RG),
                 637529648.5/(datos.acqus.NS*datos.acqus.RG)]  # tercer pico, aproximado
    center = [-0.7, 3.617, 7.347] # ppm
    fwhm = [2341.7218 / datos.acqus.SFO1,
            1873.5914 / datos.acqus.SFO1,
            2663.8478 / datos.acqus.SFO1]  # tercer pico
    fraction = [0.1379, 0.5785, 0.6014]  # Lorentzian fraction

    # Ajuste pseudo-Voigt
    vfit = PseudoVoigtFit(xdata,
                          ydata,
                          Npicos=Npicos,
                          ajustar=True,
                          amplitude=amplitude,
                          center=center,
                          fwhm=fwhm,
                          fraction=fraction,
                          fijar=["center", "fwhm", "fraction"]
                          )

    # Extraer amplitudes, errores y centros
    areas = [vfit.params[f"m{i+1}_amplitude"].value for i in range(Npicos)]
    areas_err = [vfit.params[f"m{i+1}_amplitude"].stderr for i in range(Npicos)]
    centers = [vfit.params[f"m{i+1}_center"].value for i in range(Npicos)]
    centers_err = [vfit.params[f"m{i+1}_center"].stderr for i in range(Npicos)]

    # Guardar resultados en la lista
    entry = {"expn": expn, "contact_time": contact_time}
    for i in range(Npicos):
        entry[f"m{i+1}_area"] = areas[i]
        entry[f"m{i+1}_area_err"] = areas_err[i]
        entry[f"m{i+1}_center"] = centers[i]
        entry[f"m{i+1}_center_err"] = centers_err[i]

    results.append(entry)

    # Mostrar ajuste individual
    if show_individual_fits:
        fig_fit = vfit.plot_ajuste(xlabel=r"$^1$H Chemical Shift [ppm]",
                                   ylabel="Intensity (a.u.)",
                                   reverse_xaxis=True)
        fig_fit.gca().set_title(f"Exp {expn} — Contact time = {contact_time:.1f} μs")

#=====================================================================
# Crear DataFrame y ordenarlo
#=====================================================================

df_results = pd.DataFrame(results)
df_results.sort_values(by="contact_time", inplace=True, ignore_index=True)
print(df_results)

#=====================================================================
# Gráfico de áreas
#=====================================================================

colors = ['r', 'b', 'g']
fig_area, ax_area = plt.subplots(num=382910)
for i in range(Npicos):
    ax_area.errorbar(df_results["contact_time"],
                     df_results[f"m{i+1}_area"],
                     yerr=1.96*df_results[f"m{i+1}_area_err"],
                     fmt='o-',
                     label=f'Peak: {center[i]}',
                     color=colors[i])
ax_area.set_xlabel("Contact time")
ax_area.set_ylabel("Amplitude (a.u.)")
ax_area.legend()
ax_area.grid(True)

#=====================================================================
# Gráfico de centros
#=====================================================================

fig_center, ax_center = plt.subplots(num=382911)
for i in range(Npicos):
    ax_center.plot(df_results["contact_time"],
                   df_results[f"m{i+1}_center"],
                   'o-',
                   label=f'Peak: {center[i]}',
                   color=colors[i])
ax_center.set_xlabel("Contact time")
ax_center.set_ylabel("Chemical Shift (ppm)")
ax_center.legend()
ax_center.grid(True)

#=====================================================================
# Guardado de resultados
#=====================================================================

if save:
    fig_area.savefig(f"{savepath}/{muestra}_Areas_vs_contactTime.png")
    fig_center.savefig(f"{savepath}/{muestra}_Centers_vs_contactTime.png")
    df_results.to_csv(f"{savepath}/{muestra}_VoigtFit_results.csv", index=False)
