# -*- coding: utf-8 -*-
"""
Adapted on Sat Nov 8 2025

Reemplaza la integración por ajuste de dos picos pseudo-Voigt.
Guarda los resultados en un DataFrame ordenado por contact time.
Adicional: guarda espectros ajustados y componentes, grafica stack.
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

# expns = np.concatenate([np.arange(60, 65), np.arange(66, 71)])
expns = [67,71,62,68,63,70,60,64,66,69]
# expns = [69]  # para prueba
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
# guess inicial
amplitude = [466021696.3/(datos.acqus.NS*datos.acqus.RG),
            604812667/(datos.acqus.NS*datos.acqus.RG)]
center = [3.6, -1.2] # ppm
fwhm = [2602.5823 / datos.acqus.SFO1,
        3220.4293 / datos.acqus.SFO1]
fraction = [0.6742, 0.0772]  # Lorentzian fraction



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


    vfit = PseudoVoigtFit(xdata,
                    ydata,
                    Npicos=2,
                    ajustar=True,
                    amplitude=amplitude,
                    center=center,
                    fwhm=fwhm,
                    fraction=fraction,
                    bounds={"m1_center": (center[0]-0.2, center[0]+0.2),
                            "m2_center": (center[1]-0.2, center[1]+0.2),
                            "m1_fraction": (fraction[0]*0.8, np.min([fraction[0]*1.2, 1])),
                            "m2_fraction": (fraction[1]*0.8, np.min([fraction[1]*1.2, 1])),
                            "m1_ffwhm": (fwhm[0]*0.8, fwhm[0]*1.2),
                            "m2_ffwhm": (fwhm[1]*0.8, fwhm[1]*1.2)
                            }
                    # fijar=["fwhm"]
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

    amplitude = [m1_area, m2_area]
    center = [m1_center, m2_center]
    fraction = [vfit.params["m1_fraction"].value, vfit.params["m2_fraction"].value]
    fwhm = [vfit.params["m1_fwhm"].value, vfit.params["m2_fwhm"].value]

    # Guardar componentes
    yfit, componentes = vfit.componentes(xdata)
    ycomp1, ycomp2 = componentes

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
        "m2_center_err": m2_center_err,
        "ydata": ydata,
        "ycomp1": ycomp1,
        "ycomp2": ycomp2,
        "yfit": yfit
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
# Gráfico de resultados: AREA
#=====================================================================
colors = ['r', 'b']
fig_area, ax_area = plt.subplots(num=382910)
ax_area.plot(df_results["contact_time"], df_results["m1_area"], 'o-', label=f'Peak: {center[0]}', color=colors[0])
ax_area.plot(df_results["contact_time"], df_results["m2_area"], 'o-', label=f'Peak: {center[1]}', color=colors[1])
ax_area.set_xlabel("Contact time")
ax_area.set_ylabel("Amplitude (a.u.)")
ax_area.legend()
ax_area.grid(True)
# ax_area.set_xscale('log')

#%%=====================================================================
# Gráfico de resultados: CENTER
#=====================================================================
fig_center, ax_center = plt.subplots(num=382911)
ax_center.plot(df_results["contact_time"], df_results["m1_center"], 'o-', label=f'Peak: {center[0]}', color=colors[0])
ax_center.plot(df_results["contact_time"], df_results["m2_center"], 'o-', label=f'Peak: {center[1]}', color=colors[1])
ax_center.set_xlabel("Contact time")    
ax_center.set_ylabel("Center (ppm)")
ax_center.legend()
ax_center.grid(True)

#%%=====================================================================
# Stack de espectros con componentes y ajuste total, mostrando contact time
#=====================================================================
fig_stack, ax_stack = plt.subplots(num=555555, figsize=(3,8))
contact_times_array = df_results["contact_time"].values
ydata_array = np.array(df_results["ydata"].to_list())
ycomp1_array = np.array(df_results["ycomp1"].to_list())
ycomp2_array = np.array(df_results["ycomp2"].to_list())
yfit_array = np.array(df_results["yfit"].to_list())

# ordenar por contact_time
sort_idx = np.argsort(contact_times_array)
contact_times_sorted = contact_times_array[sort_idx]
ydata_sorted = ydata_array[sort_idx, :]
ycomp1_sorted = ycomp1_array[sort_idx, :]
ycomp2_sorted = ycomp2_array[sort_idx, :]
yfit_sorted = yfit_array[sort_idx, :]

for i, t in enumerate(contact_times_sorted):
    offset = 1
    n = np.max(ydata_sorted[i]) 
    # Graficar espectros
    ax_stack.plot(ppmAxis[mask], ydata_sorted[i]/n + i*offset, color='grey', lw=4)
    ax_stack.plot(ppmAxis[mask], yfit_sorted[i]/n + i*offset, color='k', lw=1)
    ax_stack.plot(ppmAxis[mask], ycomp1_sorted[i]/n + i*offset, color=colors[0], lw=0.8)
    ax_stack.plot(ppmAxis[mask], ycomp2_sorted[i]/n + i*offset, color=colors[1], lw=0.8)
    
    # Agregar texto con contact time en microsegundos
    text_x = np.max(ppmAxis[mask]) - 35  # un poco a la izquierda
    ax_stack.text(text_x, i*offset + 0.7, f"{t:.0f} µs", va='center', ha='right', fontsize=12)

ax_stack.set_xlabel("Chemical Shift (ppm)")
ax_stack.set_ylabel("Stacked Spectra")
ax_stack.set_xlim(np.max(ppmAxis[mask]), np.min(ppmAxis[mask]))
ax_stack.grid(True)


#%%=====================================================================
# Gráfico de espectros superpuestos normalizados (solo datos)
# con colormap viridis y colorbar discreta con labels por contact time
#=====================================================================

import matplotlib.colors as mcolors

fig_super, ax_super = plt.subplots(figsize=(6,4), num=999999)

# Extraer datos
contact_times_array = df_results["contact_time"].values
ydata_array = np.array(df_results["ydata"].to_list())
ppmAxis_masked = ppmAxis[mask]

# Normalizar cada espectro por su máximo
ydata_norm = ydata_array / ydata_array.max(axis=1)[:, None]

# Generar una lista de colores, uno por espectro
num_spectra = len(contact_times_array)
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, num_spectra))

# Graficar espectros
offset = 0.2
for i, (t, color) in enumerate(zip(contact_times_array, colors)):
    ax_super.plot(ppmAxis_masked, ydata_norm[i]+offset*i, color=color, lw=3, alpha=0.9)
# Formato del gráfico
ax_super.set_xlim(np.max(ppmAxis_masked), np.min(ppmAxis_masked))
ax_super.set_xlabel("Chemical Shift (ppm)")
ax_super.set_ylabel("Normalized Intensity (a.u.)")
# ax_super.axvline(-1.2, color='k', lw=1, ls='--')
# ax_super.axvline(3.7, color='k', lw=1, ls='--')
ax_super.grid(True)

# Crear colorbar discreta con un color por espectro
cmap_discrete = mcolors.ListedColormap(colors)
bounds = np.arange(num_spectra + 1)
norm = mcolors.BoundaryNorm(bounds, cmap_discrete.N)
sm = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm)
sm.set_array([])

# Agregar colorbar con etiquetas correspondientes a contact times
cbar = fig_super.colorbar(sm, ax=ax_super, ticks=bounds[:-1] + 0.5)
cbar.ax.set_yticklabels([f"{int(t)}" for t in contact_times_array])
cbar.set_label("Contact time (µs)")

#=====================================================================
# Guardado de resultados
#=====================================================================
if save:
    fig_area.savefig(f"{savepath}/{muestra}_Areas_vs_contactTime.png")
    fig_center.savefig(f"{savepath}/{muestra}_Centers_vs_contactTime.png")
    fig_stack.savefig(f"{savepath}/{muestra}_StackSpectra.png")
    df_results.to_csv(f"{savepath}/{muestra}_VoigtFit_results.csv", index=False)

# %%
