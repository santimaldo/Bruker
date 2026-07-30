# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

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
                x = eval(contenido)
                return np.array(x)
    raise ValueError("No se encontro una linea que empiece con 'x='.")

def shift_spectrum_by_bsms(spectrum, ppmAxis, bsms, shift_factor=-1.1665e-2):
    """
    Desplaza un espectro en puntos enteros de la discretización ppm.

    Parámetros
    ----------
    spectrum : np.array
        Espectro 1D (real o complejo).
    ppmAxis : np.array
        Eje ppm original.
    bsms : float
        Valor de bsms correspondiente al espectro.
    shift_factor : float
        Factor de conversión bsms -> ppm.

    Retorna
    -------
    shifted_spectrum : np.array
        Espectro desplazado en puntos enteros.
    """
    # resolución digital del eje ppm
    dppm = np.abs(ppmAxis[1] - ppmAxis[0])
    # corrimiento continuo esperado
    shift_ppm = shift_factor * bsms
    # cantidad entera de puntos digitales
    n_shift = int(np.round(shift_ppm / dppm))
    # desplazamiento sin condiciones periódicas
    shifted_spectrum = np.zeros_like(spectrum)
    if n_shift > 0:
        shifted_spectrum[n_shift:] = spectrum[:-n_shift]
    elif n_shift < 0:
        shifted_spectrum[:n_shift] = spectrum[-n_shift:]
    else:
        shifted_spectrum = spectrum.copy()
    return shifted_spectrum

################## end Functions #######################


absolute = False
autoph = False
path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_08_08_2025_gf410_7/"
# directorio de guradado
savepath= path
expns = [12] # [5, 8, 9, 12]
muestra = "gf410_7"
save = False
plotRange = [300,-30]
# rango de integracion
ppmRanges = [
             [400, 15],
             [15, -15]            
            ]

absolute = False
autoph = False
path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_07_08_2025_gf410_5/"
# directorio de guradado
savepath= path
expns = [12, 15] #
muestra = "gf410_5"
save = False
plotRange = [300,-30]
# rango de integracion
ppmRanges = [[400, 15],
             [15, -15]            
            ]




# absolute = False
# autoph = False
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_16_09_2025_gf410_23/"
# # directorio de guradado
# savepath= path
# expns = [8] # [6, 8, 10]
# muestra = "gf410_23"
# save = False
# plotRange = [300,-30]
# # rango de integracion
# ppmRanges = [[400, 15],
#              [15, -15]            
#             ]


# absolute = False
# autoph =  False
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\Gabriel\Gabriel_16_09_2025_gf410_cell/"
# # directorio de guradado
# savepath= path
# expns = [29] # [27, 28, 29] ## 2.5 W
# expns = [37] # [35,37] ## 4.0 W 
# expns = [14] # [11, 12, 14] ## 5W 
# # expns = [29,37,14]
# muestra = "gf410_15?"
# save = False
# plotRange = [300,-30]
# # rango de integracion
# ppmRanges = [[400, 15],
#              [15, -15],             
#             ]

plot_together =  True
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


    Signals = np.array([])

    # grafico todos los espectros juntos
    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados2D(f'{path}/{expn}/',
                                read_pp = False)
    spec = datos.espectro.real
    speci = datos.espectro.imag
    maxx = np.max(spec)
    bsms_field = read_bsms_field(path_2D)
    norm = mpl.colors.Normalize(vmin=min(bsms_field), vmax=max(bsms_field))
    signal_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    center_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    width_perBsms = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    colors = ['k', 'b', 'r', 'forestgreen', 'cyan', 'magenta']
    integrals[expn] = dict()
    # Plot the 1D spectra
    for kk in range(len(bsms_field)):
        bsms = bsms_field[kk]
        ###ppmAxis = datos.espectro.ppmAxis - (-1.1665e-2) * bsms_field[kk] ### old. not used
        ppmAxis = datos.espectro.ppmAxis
        spec1d = shift_spectrum_by_bsms(spec[kk,:],ppmAxis,bsms_field[kk])
        speci1d = shift_spectrum_by_bsms(speci[kk,:],ppmAxis,bsms_field[kk])

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
        integrals[expn]["integrals"] = signal_perBsms
        integrals[expn]["bsms"] = bsms_field  
        integrals[expn]["center"] = center
        integrals[expn]["width"] = width
        abs_s = np.abs(signal_perBsms[:,ii])
        bsmsofmax = bsms_field[abs_s==np.max(abs_s)]
        bsmsofmincenter = bsms_field[center_perBsms[:,ii]==np.min(center_perBsms[:,ii])]
        integrals[expn]["bsmsofmax"] = bsmsofmax
        integrals[expn]["bsmsofmincenter"] = bsmsofmincenter

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
            axs_ranges[0, ii].set_title(f"bsms of:  max={bsmsofmax}, min_centre={bsmsofmincenter}")
            
    
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
bsms_max = []
for expn in expns:
    fig, ax = plt.subplots()
    fig_inv, ax_inv = plt.subplots()
    for ii,ppmRange in enumerate(ppmRanges):
        color = cmap(norm_ppmRange(ii))
        bsms = integrals[expn]["bsms"] #- integrals[expn]["bsmsofmax"]
        signal = np.real(integrals[expn]["integrals"][:,ii])
        ax.plot(bsms, signal, 'o-', color=color)
        ax_inv.plot(bsms, 1/signal, 'o-', color=color)
        bsms_max.append(bsms[signal==np.max(signal)])
    # ax.set_xlim(2000,-2000)
    # ax.set_xlabel(r"$\Delta$bsms")
    ax.set_xlabel(r"bsms")
    ax.set_ylabel("Spectrum integral in a range [a.u]")
    # colorbar usando la normalización por índice
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm_ppmRange)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    # ticks en posiciones de índice, labels en ppm
    cbar.set_ticks(range(len(ppmRanges)))
    ppm_ticks = [ppmRange[0] for ppmRange in ppmRanges]
    cbar.set_ticklabels(ppm_ticks)
    cbar.set_label("ppm")

#%%
vmax = 0.5 # Color vmax [0-1]
# Transformación bsms -> escala del segundo eje
factor = (-1.1665e-2*1e-6*9.42)*1e3
def bsms_to_secondary(y):
    return (y-bsms0) * factor
def secondary_to_bsms(y):
    return (y-bsms0) / factor
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
    spec_2D = np.array([np.real(spectra[expn][kk]["spec"]) for kk in kk_sorted])   # (N_bsms, N_ppm)
    bsms_1D = np.array([spectra[expn][kk]["bsms"] for kk in kk_sorted])        # (N_bsms,)
    bsms_2D = np.tile(bsms_1D[:, None], (1, ppm_2D.shape[1]))                  # (N_bsms, N_ppm)

    fig_map, ax_map = plt.subplots(figsize=(7, 5))
    pcm = ax_map.pcolormesh(ppm_2D, bsms_2D, spec_2D,
                            shading='auto', cmap='inferno',
                            vmax=np.max(spec_2D)*vmax)
    ax2 = ax_map.secondary_yaxis('right',functions=(bsms_to_secondary,secondary_to_bsms))
    ax2.set_ylabel(r"$\Delta B$ (mT)", labelpad=-3)
    ax_map.set_xlim(max(plotRange), min(plotRange))       
    ax_map.set_xlabel(r'$\delta$ (ppm)')
    ax_map.set_ylabel('bsms')
    ax_map.set_title(f'cell: {muestra}; expn {expn}')
    fig_map.colorbar(pcm, ax=ax_map, label='NMR signal [a.u.]', pad=0.12)

#%%================================================================
#==================================================================
#==================================================================
#==================================================================
#==================================================================
vmax = 1
factor = (-1.1665e-2*1e-6*9.42)*1e3

def bsms_to_secondary(y):
    return (y-bsms0)*factor

def secondary_to_bsms(y):
    return y/factor + bsms0

for expn in expns:
    bsms = integrals[expn]["bsms"]
    sig = np.abs(integrals[expn]["integrals"])
    # bsms0 = máximo de la primera ventana de integración
    bsms0 = bsms[np.argmax(sig[:,0])]

    kk = sorted(spectra[expn].keys())
    ppm_2D = np.array([spectra[expn][k]["ppm"] for k in kk])
    spec_2D = np.array([np.real(spectra[expn][k]["spec"]) for k in kk])
    bsms_1D = np.array([spectra[expn][k]["bsms"] for k in kk])
    bsms_2D = np.tile(bsms_1D[:,None], (1,ppm_2D.shape[1]))

    fig, (ax_int, ax_map) = plt.subplots(1,2,figsize=(10,5),
                                         sharey=True,
                                         gridspec_kw={'width_ratios':[1,3]})

    # for ii, ppmRange in enumerate(ppmRanges):
    ppmRange = ppmRanges[0]
    ax_int.plot(sig[:,ii], bsms, 'o-', color="k",#cmap(norm_ppmRange(ii)),
                label=f"{ppmRange[0]} to {ppmRange[1]} ppm")
    ax_int.set(xlabel="Integral [a.u.]", ylabel="bsms")
    ax_int.grid(alpha=.3)
    ax_int.legend(fontsize=8)

    ax_int2 = ax_int.secondary_yaxis('right',
                                     functions=(bsms_to_secondary,
                                                secondary_to_bsms))
    ax_int2.set_ylabel(r"$\Delta B$ (mT)")

    pcm = ax_map.pcolormesh(ppm_2D, bsms_2D, spec_2D,
                            shading='auto', cmap='inferno',
                            vmax=np.max(spec_2D)*vmax)
    ax_map.set(xlim=(max(plotRange),min(plotRange)),
               xlabel=r"$\delta$ (ppm)")
    ax_map.set_ylabel("bsms")
    ax_map.set_title(f'cell: {muestra}; expn {expn}')
    ax_map.tick_params(labelbottom=True, labelleft=True)
    ax_map2 = ax_map.secondary_yaxis('right',
                                     functions=(bsms_to_secondary,
                                                secondary_to_bsms))
    ax_map2.set_ylabel(r"$\Delta B$ (mT)", labelpad=-3)
    
    fig.colorbar(pcm, ax=ax_map, label="NMR signal [a.u.]", pad=.12)
    fig.tight_layout()
# %%
