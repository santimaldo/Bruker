# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

import ast  # FIX (e): reemplaza eval() por ast.literal_eval(), mas seguro
import nmrglue as ng
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from Datos import *
from Autophase import autophase
from scipy.integrate import simpson
import re

################## Functions ###########################

def read_bsms_field(path_archivo):
    """
    Lee un archivo con formato específico y devuelve el array 'x' como lista de floats.

    Parámetros:
    path_archivo (str): Ruta al archivo de entrada.

    Retorna:
    list[float]: Lista con los valores del array x.
    """
    with open(path_archivo+'Klog_opt', 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                x = ast.literal_eval(contenido)  # FIX (e): antes era eval(contenido)
                return np.array(x)
    raise ValueError("No se encontro una linea que empiece con 'x='.")

################## end Functions #######################
# FIX (d): se elimino un "return np.array(data)" muerto (referenciaba una
# variable 'data' inexistente y era inalcanzable porque el raise de arriba
# siempre corta la ejecucion antes de llegar ahi). Resto de copy/paste.


# absolute = False
# autoph = False
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_08_08_2025_gf410_7/"
# # directorio de guradado
# savepath= path
# exp_to_CalibrateTheOff_info = {"expn": 5,
#                                 "ppmRange": [400, 15], # to find the minimum
#                                 }      
# expns = [9,12]
# muestra = "gf410_7"
# save = False
# plotRange = [300,-20]
# ppm_bin_width = 10


# exp_to_CalibrateTheOff_info = {"expn": 12,
#                                 "ppmRange": [400, 15], # to find the minimum
#                                 }  
# absolute = False
# autoph = False
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_07_08_2025_gf410_5/"
# # directorio de guradado
# savepath= path
# expns = [12, 15] #
# muestra = "gf410_5"
# save = False
# plotRange = [700,-300]
# # rango de integracion
# ppmRanges = [[400, 15],
#              [15, -15]            
#             ]
# ppm_bin_width = 10

# exp_to_CalibrateTheOff_info = {"expn": 6,
#                                 "ppmRange": [400, 15], # to find the minimum
#                                 }  
# absolute = False
# autoph = False
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_16_09_2025_gf410_23/"
# # directorio de guradado
# savepath= path
# expns = [10] # [6, 8, 10] ############ tienen distintos sr!!!!!
# muestra = "gf410_23"
# save = False
# plotRange = [700,-300]
# # rango de integracion
# ppmRanges = [[400, 15],
#              [15, -15]            
#             ]
# ppm_bin_width = 10

absolute = False
autoph =  False
path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_16_09_2025_gf410_cell/"
# directorio de guradado
savepath= path
# expns = [29] # [27, 28, 29] ## 2.5 W
# expns = [37] # [35,37] ## 4.0 W 
expns = [14] # [11, 12, 14] ## 5W 
expns = [11, 12, 14]  
exp_to_CalibrateTheOff_info = {"expn": 12, # 12->5W
                                "ppmRange": [400, 15], # to find the minimum
                                }  
muestra = "gf410_15?"
save = False
plotRange = [300,-20]
ppm_bin_width = 20


plot_together =  True
#=====================================================================
# Calibrate the "off resonance"
#=====================================================================

expn = exp_to_CalibrateTheOff_info["expn"]
ppmRange = exp_to_CalibrateTheOff_info["ppmRange"]
# grafico todos los espectros juntos
path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(f'{path}/{expn}/',
                            read_pp = False)
bsms_field = read_bsms_field(path_2D)
spec = datos.espectro.real
speci = datos.espectro.imag
# Plot the 1D spectra
specs = []
integrals_ref = []
for kk in range(len(bsms_field)):    
    ppmAxis = datos.espectro.ppmAxis - (-1.1665e-2) * bsms_field[kk]
    spec1d = spec[kk,:]
    speci1d = speci[kk,:]
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
    # signal = datos.Integrar(ppmRange=ppmRange)
    condition = np.abs(ppmAxis-np.mean(ppmRange)) < np.abs(r1-r2)/2
    spec_to_integrate = (spec1d+1j*speci1d)[condition]
    ppm_to_integrate = ppmAxis[condition]
    specs.append([ppmAxis, spec1d+1j*speci1d])
    integrals_ref.append(-simpson(spec_to_integrate, x=ppm_to_integrate))
integrals_ref = np.array(integrals_ref)
idx_min = np.where(integrals_ref.real == np.min(integrals_ref.real))[0][0]
ppm, spec = specs[idx_min]
ppmofmax = ppm[spec.real == np.max(spec.real)][0]
condition = np.abs(ppm-ppmofmax) < ppm_bin_width/2 # elijo un rango de +- 10 ppm 
ppm_to_integrate =ppm[condition]
spec_to_integrate = spec[condition]
integral_off = -simpson(spec_to_integrate, x=ppm_to_integrate)
fig, ax = plt.subplots()
ax.plot(ppm, spec.real)
ax.plot(ppm, spec.imag, '--', color="tab:blue")
ax.axvspan(ppmofmax-ppm_bin_width/2, ppmofmax+ppm_bin_width/2, alpha=0.2)
ax.set_title(f"max: {ppmofmax:.2f} ppm; integral {integral_off:.2e} (a.u)")
ax.set_xlim(max(ppmRange), min(ppmRange))


off_res_params = {"ppm": ppmofmax,
                  "integral": integral_off,
                  }
ppm_bins = np.arange(ppmofmax, ppmofmax-240, -ppm_bin_width)
# FIX (a): el ancho total de cada bin tiene que ser ppm_bin_width (no
# 2*ppm_bin_width). Antes: [ppm+ppm_bin_width, ppm-ppm_bin_width], que da
# bins de 2*ppm_bin_width de ancho, espaciados solo ppm_bin_width entre si
# -> 50% de solapamiento entre bins consecutivos. Con /2 los bins quedan
# contiguos y sin superposicion, tal como en la seccion de calibracion de
# arriba (que ya usaba correctamente ppm_bin_width/2).
ppmRanges = [[ppm+ppm_bin_width/2, ppm-ppm_bin_width/2]
             for ppm in ppm_bins]
s_bins = [1-(ppm/ppmofmax) for ppm in ppm_bins]


#=====================================================================
# 2D experiments
#=====================================================================
if plot_together:
    fig_popt, axs_popt = plt.subplots(nrows=len(ppmRanges))#,figsize=(4,2*len(ppmRanges)))
    fig_ranges, axs_ranges = plt.subplots(nrows=3, ncols=len(ppmRanges))                                          
cmap = plt.get_cmap("gnuplot")
norm_ppmRange = mpl.colors.Normalize(vmin=0, vmax=len(ppmRanges)-1)
bsms_list = []
spectra = dict.fromkeys(expns)
integrals = dict.fromkeys(expns)
for expn in expns:
    spectra[expn] = dict()
    fig_spec, ax_spec = plt.subplots(figsize=(4,15))
    if not plot_together:
        fig_popt, axs_popt = plt.subplots(nrows=len(ppmRanges))
        fig_ranges, axs_ranges = plt.subplots(nrows=3, ncols=len(ppmRanges))
    ax_spec.set_xlim(max(plotRange), min(plotRange))

    # grafico todos los espectros juntos
    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados2D(f'{path}/{expn}/',
                                read_pp = False)
    spec = datos.espectro.real
    speci = datos.espectro.imag
    maxx = np.max(spec)
    bsms_field = read_bsms_field(path_2D)

    # FIX (f): el calculo de centroide mas abajo (simpson con x=bsms) solo
    # tiene sentido si bsms_field es monotono. Se chequea aca, apenas se lee.
    diffs = np.diff(bsms_field)
    assert np.all(diffs > 0) or np.all(diffs < 0), (
        f"bsms_field no es monotono para expn={expn}; "
        "el calculo de centroide (means) mas abajo va a dar resultados "
        "sin sentido si no se ordena antes."
    )

    norm = mpl.colors.Normalize(vmin=min(bsms_field), vmax=max(bsms_field))
    signal_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    center_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    width_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
    integrals[expn] = dict()
    # Plot the 1D spectra
    for kk in range(len(bsms_field)):
        bsms = bsms_field[kk]
        ppmAxis = datos.espectro.ppmAxis - (-1.1665e-2) * bsms_field[kk]
        spec1d = spec[kk,:]
        speci1d = speci[kk,:]
        spectra[expn][kk] = dict()
        spectra[expn][kk]["spec"] = spec1d + 1j* speci1d 
        spectra[expn][kk]["ppm"] = ppmAxis
        spectra[expn][kk]["bsms"] = bsms
        
        if absolute:
            abs_spec = np.abs(spec1d + 1j*speci1d)
            spec[kk,:]  = abs_spec
            spec1d = abs_spec
        elif autoph:
            spec1d = ng.proc_autophase.autops(spec1d+1j*speci1d, "acme") #  "peak_minima")
            spec[kk,:]  = ng.process.proc_bl.cbf(spec1d.real, last=100)
        ax_spec.axhline(maxx * kk, ls='--', color='grey', lw=2) 
        ax_spec.plot(ppmAxis, spec1d+maxx * kk, color = 'k')# cmap(norm(bsms)))
        ax_spec.plot(ppmAxis, speci1d+maxx * kk, '--', color = 'k')#cmap(norm(bsms)))
        ax_spec.text(max(plotRange)*0.95, maxx*(kk+0.2), bsms)#f"bsms: {bsms}")
            
        ###### start integrating the spectra=
        for ii, ppmRange in enumerate(ppmRanges):
            color = cmap(norm_ppmRange(ii))
            r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
            if kk==len(bsms_field)-1: # solo grafico pal ultimo
                ax_spec.axvspan(r1, r2, alpha=0.1, color=color)

            # signal = datos.Integrar(ppmRange=ppmRange)
            condition = np.abs(ppmAxis-np.mean(ppmRange)) < np.abs(r1-r2)/2
            spec_to_integrate = (spec1d+1j*speci1d)[condition]
            ppm_to_integrate = ppmAxis[condition]  
            # los doy vuelta para calcular los momentos de manera
            # mas prolija
            ppm_to_integrate = ppm_to_integrate[::-1]
            spec_to_integrate = spec_to_integrate[::-1]   
            # Integral (área)
            y = spec_to_integrate
            x = ppm_to_integrate
            signal_perBsms[kk, ii] = simpson(y, x=x)
            # Centro (primer momento)
            y = np.abs(spec_to_integrate)
            center = simpson(ppm_to_integrate * y,x=x)
            center/= simpson(y, x=x)
            # Ancho RMS (segundo momento central)
            width = np.sqrt(
                simpson((ppm_to_integrate-center)**2 * y, x=x) / 
                simpson(y, x=x))
            center_perBsms[kk, ii] = center
            width_perBsms[kk, ii] = width

    # FIX (b) y (c): este bloque estaba antes DENTRO del "for kk" (se
    # recalculaba, y se pisaba, en cada iteracion de kk -- no rompia el
    # resultado final pero era trabajo redundante). Ahora corre UNA vez,
    # despues de que signal_perBsms/center_perBsms/width_perBsms ya estan
    # completamente llenos para todos los kk.
    #
    # Ademas, antes "ii" quedaba pegado al ULTIMO valor del loop interno
    # "for ii, ppmRange in enumerate(ppmRanges)", entonces bsmsofmax /
    # bsmsofmincenter (y "center"/"width") terminaban calculados solo para
    # el ULTIMO bin de ppm, no para todos. Ahora se calculan por bin
    # (arrays de largo len(ppmRanges)), y se guardan los arrays completos
    # center_perBsms / width_perBsms en vez de un escalar suelto.
    integrals[expn]["integrals"] = signal_perBsms
    integrals[expn]["bsms"] = bsms_field
    integrals[expn]["center_perBsms"] = center_perBsms
    integrals[expn]["width_perBsms"] = width_perBsms

    bsmsofmax_perBin = np.zeros(len(ppmRanges))
    bsmsofmincenter_perBin = np.zeros(len(ppmRanges))
    for ii in range(len(ppmRanges)):
        abs_s = np.abs(signal_perBsms[:, ii])
        bsmsofmax_perBin[ii] = bsms_field[abs_s == np.max(abs_s)][0]
        bsmsofmincenter_perBin[ii] = bsms_field[
            center_perBsms[:, ii].real == np.min(center_perBsms[:, ii].real)
        ][0]
    integrals[expn]["bsmsofmax"] = bsmsofmax_perBin
    integrals[expn]["bsmsofmincenter"] = bsmsofmincenter_perBin

    for ii in range(len(ppmRanges)):
        color = cmap(norm_ppmRange(ii))
        ax_popt = axs_popt[ii]
        ax_popt.plot(bsms_field, np.real(signal_perBsms[:,ii]), 'bo-')
        ax_popt.plot(bsms_field, np.imag(signal_perBsms[:,ii]), 'ro-')
        ax_popt.plot(bsms_field, np.abs(signal_perBsms[:,ii]), 'ko-')
        
        
        axs_ranges[0, ii].plot(bsms_field, np.abs(signal_perBsms[:,ii]), 'o-', color=color)
        axs_ranges[1, ii].plot(bsms_field, center_perBsms[:,ii], 'o-', color=color)
        axs_ranges[2, ii].plot(bsms_field, width_perBsms[:,ii], 'o-', color=color)

        if ii == 0:
            # FIX (b): ahora si corresponde al bin ii=0 (antes mostraba
            # los valores del ultimo bin, quedaban de un "ii" viejo)
            axs_ranges[0, ii].set_title(
                f"bsms of:  max={bsmsofmax_perBin[ii]:.1f}, "
                f"min_centre={bsmsofmincenter_perBin[ii]:.1f}"
            )
            
    
#     if save:        
#         for kk in range(len(bsms_field)):
#             np.savetxt(#f"{savepath}/{muestra}_expn{expn}_bsms_{int(bsms_field[kk]):04d}.dat",
#                 f"{savepath}/{muestra}_bsms_{int(bsms_field[kk]):04d}.dat",
#                     np.array([ppmAxis, spec[kk,:]]).T,
#                     header="ppmAxis\treal")
#     # ax_popt.plot(bsms_field, Signals, 'o')#, color=color, label="Rising Edge")
# np.savetxt(f"{savepath}/bsms_list",
#            bsms_list)
# print(np.array([all_bsms_field, all_Signals]).T)

# %%
# lists_of_plots = []
# # lists_of_plots.append([5, 2]) # expn 29
# # lists_of_plots.append([6, 2]) # expn 37
# lists_of_plots.append([28,19]) # expn 14
# # fig, ax = plt.subplots()
# for ii,expn in enumerate(expns):
#     spec = spectra[expn]
#     fig, ax = plt.subplots()
#     for idx in lists_of_plots[ii]:
#         line, =ax.plot(spec[idx]["ppm"], spec[idx]["spec"].real, label=spec[idx]["bsms"])
#         ax.plot(spec[idx]["ppm"], spec[idx]["spec"].imag, '--', color=line.get_color())
#         ax.plot(spec[idx]["ppm"], np.abs(spec[idx]["spec"]), color=line.get_color(), alpha=0.3)
#     ax.legend()
#     ax.set_xlim(400, -50)
#     # ax.set_ylim([-0.2e7, 6e7])



# %%
ax_dens_ylim_upper = 0.15e6
integral_off = off_res_params["integral"]                  
emax = 1694
for expn in expns:
    fig, ax = plt.subplots()
    fig_dens, ax_dens = plt.subplots()
    means = []
    maxss = []
    bsms = integrals[expn]["bsms"] #- integrals[expn]["bsmsofmax"]
    for ii,ppmRange in enumerate(ppmRanges):
        color = cmap(norm_ppmRange(ii))        
        signal = np.real(integrals[expn]["integrals"][:,ii])
        signal_off = 1 # np.real(integral_off)
        print("WARNING! ESTOY USANDO signal_off=1")
        ax.plot(bsms, signal, 'o-', color=color)        
        density = signal/signal_off/(1+emax*s_bins[ii])
        ax_dens.plot(bsms, density, 'o-', color=color)                
        means.append(simpson(density*bsms, x=bsms)/simpson(density, x=bsms))
        maxss.append(bsms[signal==np.max(signal)])
    # ax.set_xlim(2000,-2000)
    # ax.set_xlabel(r"$\Delta$bsms")
    ax.set_xlabel(r"bsms")
    ax.set_ylabel("Spectrum integral in a range [a.u]")
    # colorbar usando la normalización por índice
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm_ppmRange)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.08)
    # ticks en posiciones de índice, labels en ppm
    cbar.set_ticks(range(len(ppmRanges)))
    ppm_ticks = [f"{ppmRange[0]:.0f}" for ppmRange in ppmRanges]
    cbar.set_ticklabels(ppm_ticks)
    cbar.set_label("ppm")
    ax2 = ax.twinx()
    signal_perBsms = np.sum(integrals[expn]["integrals"], axis=1)
    color_ax2="red"
    ax2.plot(bsms, np.real(signal_perBsms), 's-', color=color_ax2, alpha=0.5)
    ax2.spines["right"].set_color(color_ax2)
    ax2.tick_params(axis="y", colors=color_ax2)    
    ax2.set_ylabel("Total integral", color=color_ax2)
    ax_dens.set_ylim([0,None])
    ax2.set_ylim([0,None])  
    #============================================================
    ax_dens.set_xlabel(r"bsms")
    ax_dens.set_ylabel("Density")
    # colorbar usando la normalización por índice
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm_ppmRange)
    sm.set_array([])    
    cbar = fig_dens.colorbar(sm, ax=ax_dens, pad=0.08)
    # ticks en posiciones de índice, labels en ppm
    cbar.set_ticks(range(len(ppmRanges)))    
    cbar.set_ticklabels([f"{v:.2f}" for v in s_bins])
    cbar.set_label("saturation factor")    
    ax2 = ax_dens.twinx()
    color_ax2="red"
    ax2.plot(bsms, np.real(signal_perBsms), 's-', color=color_ax2, alpha=0.5)
    ax2.spines["right"].set_color(color_ax2)
    ax2.tick_params(axis="y", colors=color_ax2)    
    ax2.set_ylabel("Total integral", color=color_ax2)
    ax_dens.set_ylim([0,ax_dens_ylim_upper])
    ax2.set_ylim([0,None])
    #=============================================================
    fig, ax = plt.subplots()
    ax.plot(s_bins, means, 'bo-')
    ax.set_ylim([min(bsms), max(bsms)])    
    # ax.plot(s_bins, maxss, 'ro-')
    
#%% densidad de spines con saturacion s:

for expn in expns:
    # ============================================================
    # Colormap: ppmAxis (eje x) vs bsms (eje y)
    # ============================================================
    # Cada fila (cada bsms) tiene su PROPIO eje de ppm, ya corregido
    # (ppmAxis = ppmAxis_crudo - (-1.1665e-2)*bsms). Para no perder esa
    # correccion ni interpolar de mas, se arma un pcolormesh con X e Y
    # como grillas 2D: cada fila de X lleva el ppmAxis que corresponde a
    # esa fila de Y (bsms), en vez de forzar un unico eje x compartido.
    kk_sorted = sorted(spectra[expn].keys())
    ppm_2D = np.array([spectra[expn][kk]["ppm"] for kk in kk_sorted])          # (N_bsms, N_ppm)
    spec_2D = np.array([spectra[expn][kk]["spec"].real for kk in kk_sorted])   # (N_bsms, N_ppm)
    bsms_1D = np.array([spectra[expn][kk]["bsms"] for kk in kk_sorted])        # (N_bsms,)
    bsms_2D = np.tile(bsms_1D[:, None], (1, ppm_2D.shape[1]))                  # (N_bsms, N_ppm)

    fig_map, ax_map = plt.subplots(figsize=(6, 5))
    pcm = ax_map.pcolormesh(ppm_2D, bsms_2D, spec_2D,
                            shading='auto', cmap='inferno',
                            vmax=np.max(spec_2D)*1)
    ax_map.set_xlim(max(plotRange), min(plotRange))   # convencion ppm descendente
    ax_map.set_xlabel(r'$\delta$ (ppm)')
    ax_map.set_ylabel('bsms')
    ax_map.set_title(f'expn {expn}')
    fig_map.colorbar(pcm, ax=ax_map, label='Señal (parte real) [a.u.]')
 # %%
