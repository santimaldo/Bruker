# -*- coding: utf-8 -*-
"""
Adapted on Sat Nov 8 2025

Reemplaza la integración por ajuste de tres picos pseudo-Voigt.
Guarda los resultados en un DataFrame ordenado por contact time.
Adicional: guarda espectros ajustados y grafica stack.
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

    # Guess inicial para los tres picos
    amplitude = [359451094.1/(datos.acqus.NS*datos.acqus.RG),
                 162708124.8/(datos.acqus.NS*datos.acqus.RG),
                 637529648.5/(datos.acqus.NS*datos.acqus.RG)]
    center = [-0.7, 3.617, 7.347] # ppm
    fwhm = [2341.7218 / datos.acqus.SFO1,
            1873.5914 / datos.acqus.SFO1,
            2663.8478 / datos.acqus.SFO1]
    fraction = [0.1379, 0.5785, 0.6014]

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

    # Guardar espectros y componentes usando componentess()
    yfit, componentes = vfit.componentes(xdata)
    ycomp1, ycomp2, ycomp3 = componentes

    # Guardar resultados en la lista
    entry = {
        "expn": expn,
        "contact_time": contact_time,
        "ydata": ydata,
        "ycomp1": ycomp1,
        "ycomp2": ycomp2,
        "ycomp3": ycomp3,
        "yfit": yfit
    }

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

#%%=====================================================================
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
ax_area.set_xscale('log')

#%%=====================================================================
# Stack de espectros con componentes y ajuste total
#=====================================================================

fig_stack, ax_stack = plt.subplots(num=555555, figsize=(3,8))
contact_times_array = df_results["contact_time"].values
ydata_array = np.array(df_results["ydata"].to_list())
ycomp1_array = np.array(df_results["ycomp1"].to_list())
ycomp2_array = np.array(df_results["ycomp2"].to_list())
ycomp3_array = np.array(df_results["ycomp3"].to_list())
yfit_array = np.array(df_results["yfit"].to_list())

# Ordenar por contact_time
sort_idx = np.argsort(contact_times_array)
contact_times_sorted = contact_times_array[sort_idx]
ydata_sorted = ydata_array[sort_idx, :]
ycomp1_sorted = ycomp1_array[sort_idx, :]
ycomp2_sorted = ycomp2_array[sort_idx, :]
ycomp3_sorted = ycomp3_array[sort_idx, :]
yfit_sorted = yfit_array[sort_idx, :]

offset = 1
for i, t in enumerate(contact_times_sorted):
    n = np.max(ydata_sorted[i])
    ax_stack.plot(ppmAxis[mask], ydata_sorted[i]/n + i*offset, color='grey', lw=4)
    ax_stack.plot(ppmAxis[mask], yfit_sorted[i]/n + i*offset, color='k', lw=1)
    ax_stack.plot(ppmAxis[mask], ycomp1_sorted[i]/n + i*offset, color=colors[0], lw=0.8)
    ax_stack.plot(ppmAxis[mask], ycomp2_sorted[i]/n + i*offset, color=colors[1], lw=0.8)
    ax_stack.plot(ppmAxis[mask], ycomp3_sorted[i]/n + i*offset, color=colors[2], lw=0.8)

ax_stack.set_xlabel("Chemical Shift (ppm)")
ax_stack.set_ylabel("Stacked Spectra")
ax_stack.set_xlim(np.max(ppmAxis[mask]), np.min(ppmAxis[mask]))
ax_stack.grid(True)

#=====================================================================
# Guardado de resultados
#=====================================================================

if save:
    fig_area.savefig(f"{savepath}/{muestra}_Areas_vs_contactTime.png")
    fig_stack.savefig(f"{savepath}/{muestra}_StackSpectra.png")
    df_results.to_csv(f"{savepath}/{muestra}_VoigtFit_results.csv", index=False)
