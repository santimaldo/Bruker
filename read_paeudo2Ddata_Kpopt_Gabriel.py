# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi

"""

# # Clear variables from a previous run when this script is executed again.
# for _name in list(globals()):
#     if not _name.startswith("__"):
#         del globals()[_name]

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import os
from Datos import *
from Autophase import autophase
from scipy.integrate import simpson
import re


################## Functions ###########################
# conversion factor:
mT_per_bsms = (-1.1665e-2*1e-6*9.42)*1e3

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
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\2026-08-17_1.3mm_Gabriel_gf410_41/"
# # sample label, expn, "0 ppm"@bsms=0 (para calibrar)
# infos = []
# infos.append(["cell41_5.2W", 33, 0])


# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\2026-08-18_1.3mm_Gabriel_gf410_38/"
# # sample label, expn, "0 ppm"@bsms=0 (para calibrar)
# infos = []
# infos.append(["cell38_5.2W", 16, 0])

# # path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\2026-08-26_1.3mm_Gabriel_gf410_46/"
# # # sample label, expn, "0 ppm"@bsms=0 (para calibrar)
# # infos = []
# # infos.append(["cell46_5.2W", 17, 105.15])
# # infos.append(["cell46_4.0W", 78, 196.16])
# # infos.append(["cell46_3.0W", 76, 196.16])
# # infos.append(["cell46_2.0W", 84, 196.16])
# # infos.append(["cell46_1.0W", 88, 196.16])

# # directorio de guradado
# savepath_base= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Gabriel\Analysis\2026-08_staticDNP/"
# plotRange = [300,-30]
# # rango de integracion
# ppmRanges = [[300, 15],
#              [300, 210],
#              [210, 215],
#              [15, -15],             
#             ]





path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\2026-09-10_3.2mm_LiMetal_LiF_RC1/"

# # directorio de guradado
# savepath_base= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\analysis\2026-09_Limetal_LiF\10kHz/"
# sample label, expn, "0 ppm"@bsms=0 (para calibrar)
infos = []
# infos.append(["kpopt_p1-4.35us_p90", 28, 0, -4750])
# infos.append(["kpopt_p1-1us", 212, 0, -4750])
# infos.append(["kpopt_p1-0.5us", 35, 0, -4750])
# fields_to_plot_mT = [0, 0.1, -0.1, -0.175, -0.3]

# directorio de guradado
savepath_base= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\analysis\2026-09_Limetal_LiF\0kHz/"
# sample label, expn, "0 ppm"@bsms=0 (para calibrar), bsms@"0 mT"
fields_to_plot_mT = [0.027, 0, -0.027, -0.083, -0.11, -0.137, -0.25, 0.15]
infos = []
infos.append(["kpopt_p1-4.35us_p90", 265, 0, -4750])
# infos.append(["kpopt_p1-1us", 262, -98, -4750]) # -3750

# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp\2026-09-16_InSitu_cellTest-DNP/"
# # directorio de guradado
# savepath_base= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\analysis\2026-09_Limetal_LiF\0kHz/"
# # sample label, expn, "0 ppm"@bsms=0 (para calibrar)
# fields_to_plot_mT = [0, -0.33]
# infos = []
# infos.append(["kpopt_expn_52", 52, 0])

plotRange = [300,-50]
# rango de integracion
ppmRanges = [[300, 10],
             [0, -50],
             [300, 230],
             [230, 10],             

            ]


# ############ Gabrie g410_7
# path  =r"C:\Users\Santi\Documents\NMRdata\400dnp/Gabriel/Gabriel_08_08_2025_gf410_7/"
# # # sample label, expn, "0 ppm"@bsms=0 (para calibrar)
# infos = []
# infos.append(["g410_7", 12, -4.61])
# # directorio de guradado
# savepath_base= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Gabriel\Analysis\2025_results/"
# plotRange = [300,-30]
# # rango de integracion
# ppmRanges = [[300, 30],
#              [30, -30],             
#             ]



for info in infos:

    muestra, expn, ppmzero_at_central_bsms, set0mT_tothis_bsms = info
    savepath = os.path.join(savepath_base, muestra)
    os.makedirs(savepath, exist_ok=True)
    savepath_spec = os.path.join(savepath, "spectra")
    os.makedirs(savepath_spec, exist_ok=True)
    #=====================================================================
    # 2D experiments
    #=====================================================================

    fig_popt, axs_popt = plt.subplots(nrows=len(ppmRanges))#,figsize=(4,2*len(ppmRanges)))
    fig_ranges, axs_ranges = plt.subplots(nrows=3, ncols=len(ppmRanges))                                          
    cmap = plt.get_cmap("gnuplot")
    norm_ppmRange = mpl.colors.Normalize(vmin=0, vmax=len(ppmRanges)-1)

    spectra = dict()
    Signals = np.array([])
    integrals = dict()
    integrals_results = dict()
    #---------IMPORT DATA--------------------
    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados2D(f'{path}/{expn}/',
                                read_pp = False)
    datos.espectro.ppmAxis -= ppmzero_at_central_bsms
    spec = datos.espectro.real
    speci = datos.espectro.imag

    bsms_field = read_bsms_field(path_2D)
    signal_BSMS = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    center_BSMS = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    width_BSMS = np.zeros((len(bsms_field), len(ppmRanges))).astype("complex")
    height_BSMS = np.zeros(len(bsms_field))
    #----------------------------------------

    for kk in range(len(bsms_field)):
        bsms = bsms_field[kk]
        ppmAxis = datos.espectro.ppmAxis
        spec1d = shift_spectrum_by_bsms(spec[kk,:],ppmAxis,bsms_field[kk])
        speci1d = shift_spectrum_by_bsms(speci[kk,:],ppmAxis,bsms_field[kk])

        spectra[kk] = dict()
        spectra[kk]["spec"] = spec1d + 1j* speci1d 
        spectra[kk]["ppm"] = ppmAxis
        spectra[kk]["bsms"] = bsms
        
        height_BSMS[kk] = np.max(spec1d)
        ###### start integrating the spectra=
        for ii, ppmRange in enumerate(ppmRanges):
            color = cmap(norm_ppmRange(ii))
            r1, r2 = [np.min(ppmRange), np.max(ppmRange)]  # redefino el rango
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
            signal_BSMS[kk, ii] = simpson(y, x=x)
            # Centro (primer momento)
            y = np.real(spec_to_integrate)
            center = simpson(ppm_to_integrate * y,x=x)
            center/= simpson(y, x=x)
            # Ancho RMS (segundo momento central)
            width = np.sqrt(
                simpson((ppm_to_integrate-center)**2 * y, x=x) / 
                simpson(y, x=x))
            center_BSMS[kk, ii] = center
            width_BSMS[kk, ii] = width
            
    integrals["integrals"] = signal_BSMS
    integrals["bsms"] = bsms_field  
    integrals["center"] = center_BSMS[:,0]
    integrals["width"] = width_BSMS[:,0]
    s = np.real(signal_BSMS[:,0])
    bsmsofmax = bsms_field[s==np.max(s)][0]
    bsmsofmincenter = bsms_field[center_BSMS[:,ii]==np.min(center_BSMS[:,ii])][0]
    bsmsofmaxheight = bsms_field[height_BSMS==np.max(height_BSMS)][0]
    integrals_results["bsmsofmax"] = bsmsofmax
    integrals_results["bsmsofmincenter"] = bsmsofmincenter
    integrals_results["bsmsofmaxheight"] = bsmsofmaxheight
    integrals_results["mT_per_bsms"] = mT_per_bsms
    for ii in range(len(ppmRanges)):
        color = cmap(norm_ppmRange(ii))
        ax_popt = axs_popt[ii]
        ax_popt.plot(bsms_field, np.real(signal_BSMS[:,ii]), 'bo-')
        ax_popt.plot(bsms_field, np.imag(signal_BSMS[:,ii]), 'ro-')
        ax_popt.plot(bsms_field, np.abs(signal_BSMS[:,ii]), 'ko-')
        
        
        axs_ranges[0, ii].plot(bsms_field, np.abs(signal_BSMS[:,ii]), 'o-', color=color)
        axs_ranges[1, ii].plot(bsms_field, center_BSMS[:,ii], 'o-', color=color)
        axs_ranges[2, ii].plot(bsms_field, width_BSMS[:,ii], 'o-', color=color)

        if ii == 0:
            axs_ranges[0, ii].set_title(f"bsms of:  max={bsmsofmax}, min_centre={bsmsofmincenter}")

    if set0mT_tothis_bsms is None:
        central_bsms = integrals_results["bsmsofmaxheight"][0]
    else:                
        central_bsms = set0mT_tothis_bsms
    integrals["bsms_centred_area"] = integrals["bsms"] - integrals_results["bsmsofmax"]
    integrals["bsms_centred_area_mT"] = integrals["bsms_centred_area"] * mT_per_bsms
    integrals["bsms_centred"] = integrals["bsms"] - central_bsms
    integrals["bsms_centred_mT"] = integrals["bsms_centred"] * mT_per_bsms
    
    #-----------------------SAVE-------------------------------------
    # Each row contains the integrations for all ppm ranges.  Keep those rows as
    # objects so pandas sees this as one 1-D column rather than a 2-D column.
    df = integrals.copy()
    df["integrals"] = np.real(integrals["integrals"][:,0])
    df = pd.DataFrame(df)
    integrals_results = pd.DataFrame([integrals_results])

    df.to_csv(os.path.join(savepath, f"integrals.csv"), index=False)
    integrals_results.to_csv(os.path.join(savepath,"integrals_results.csv"), index=False)

    for kk, spec in enumerate(spectra.items()):
        data = np.array([spectra[kk]['ppm'],
                        np.real(spectra[kk]['spec']),
                        np.imag(spectra[kk]['spec'])]).T
        header = f"# BSMS: {df.bsms_centred[kk]}, DeltaB: {df.bsms_centred_mT[kk]}"
        filename = os.path.join(savepath_spec, f"spec{kk:03d}.dat")
        np.savetxt(filename, data, header=header)
    #----------------------------------------------------------------

    fig, ax = plt.subplots()
    fig_inv, ax_inv = plt.subplots()
    for ii,ppmRange in enumerate(ppmRanges):
        color = cmap(norm_ppmRange(ii))
        bsms = integrals["bsms"] #- integrals["bsmsofmax"]
        signal = np.real(integrals["integrals"][:,ii])
        ax.plot(bsms, signal, 'o-', color=color)
        ax_inv.plot(bsms, signal/max(signal), 'o-', color=color)
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

    vmax = 0.5 # Color vmax [0-1]
    # Transformación bsms -> escala del segundo eje
    def bsms_to_secondary(y):
        return (y-central_bsms) * mT_per_bsms
    def secondary_to_bsms(y):
        return (y-central_bsms) / mT_per_bsms
    # ============================================================
    # Colormap: ppmAxis (eje x) vs bsms (eje y)
    # ============================================================
    # Cada fila (cada bsms) tiene su PROPIO eje de ppm, ya corregido
    # (ppmAxis = ppmAxis_crudo - (-1.1665e-2)*bsms). Para no perder esa
    # correccion ni interpolar de mas, se arma un pcolormesh con X e Y
    # como grillas 2D: cada fila de X lleva el ppmAxis que corresponde a
    # esa fila de Y (bsms), en vez de forzar un unico eje x compartido.
    kk_sorted = sorted(spectra.keys())
    ppm_2D = np.array([spectra[kk]["ppm"] for kk in kk_sorted])          # (N_bsms, N_ppm)
    spec_2D = np.array([np.real(spectra[kk]["spec"]) for kk in kk_sorted])   # (N_bsms, N_ppm)
    bsms_1D = np.array([spectra[kk]["bsms"] for kk in kk_sorted])        # (N_bsms,)
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
   
    bsms = integrals["bsms"]
    sig = np.abs(integrals["integrals"])

    kk = sorted(spectra.keys())
    ppm_2D = np.array([spectra[k]["ppm"] for k in kk])
    spec_2D = np.array([np.real(spectra[k]["spec"]) for k in kk])
    bsms_1D = df["bsms_centred_mT"].values
    bsms_2D = np.tile(bsms_1D[:,None], (1,ppm_2D.shape[1]))

    # central_bsms = central_bsms
    # def bsms_to_secondary(y):
    #     return (y-central_bsms)*factor
    # def secondary_to_bsms(y):
    #     return y/factor + central_bsms


    fig, (ax_int, ax_map) = plt.subplots(1,2,figsize=(10,5),
                                            sharey=True,
                                            gridspec_kw={'width_ratios':[1,3]})


    ax_int.plot(df['integrals'], df["bsms_centred_mT"], 'o-', color="k",#cmap(norm_ppmRange(ii)),
               label=f"integral {ppmRanges[0][0]} to {ppmRanges[0][1]} ppm")
    print("Aca esta el plot temporal")
    # ax_int.plot(20*integrals['integrals'][:,1], df["bsms_centred_mT"], 'o-', color="k", markerfacecolor="white",#cmap(norm_ppmRange(ii)),
    #                    label=fr"20$\times$integral {ppmRanges[1][0]} to {ppmRanges[1][1]} ppm")
    ax_int.plot(integrals['integrals'][:,2], df["bsms_centred_mT"], '^-', color="k",#cmap(norm_ppmRange(ii)),
                   label=f"integral {ppmRanges[2][0]} to {ppmRanges[2][1]} ppm")
    ax_int.plot(integrals['integrals'][:,3], df["bsms_centred_mT"], 'x-', color="k",#cmap(norm_ppmRange(ii)),
                       label=f"integral {ppmRanges[3][0]} to {ppmRanges[3][1]} ppm")
    ax_int.set(xlabel="Integral [a.u.]", ylabel=r"$\Delta B$ (mT)")
    # for kk in range(len(ppmRanges)):
    #     integs = np.real(integrals['integrals'][:,kk])
    #     ax_int.plot(integs/max(integs), df["bsms_centred_mT"], 'o-',#cmap(norm_ppmRange(ii)),
    #                 label=f"integral {ppmRanges[kk][0]} to {ppmRanges[kk][1]} ppm")
    # ax_int.set(xlabel="normalized integral", ylabel=r"$\Delta B$ (mT)")
    
    ax_int.grid(alpha=.3)
    ax_int.legend(fontsize=8, loc="upper right")

    pcm = ax_map.pcolormesh(ppm_2D, bsms_2D, spec_2D,
                            shading='auto', cmap='inferno',
                            vmax=np.max(spec_2D)*vmax)
    ax_map.set(xlim=(max(plotRange),min(plotRange)),
            ylim=(max(bsms_1D),min(bsms_1D)),
            xlabel=r"$\delta$ (ppm)")
    ax_map.set_ylabel(r"$\Delta B$ (mT)")
    ax_map.set_title(f'cell: {muestra}; expn {expn}')
    ax_map.tick_params(labelbottom=True, labelleft=True)
    ax_map.axhline(0, color="white", alpha=0.75, ls="--", lw=0.5)

    fig.colorbar(pcm, ax=ax_map, label="NMR signal [a.u.]")


    fig.tight_layout()
    fig.savefig(os.path.join(savepath, "2dMap_figure.png"))

    ax_map.set_ylim(0.2, -0.4) # este solia ser hasta -0.2
    ax_int.set_ylim(0.2, -0.4) # este solia ser hasta -0.2
    fig_int = fig
    

    # espectros a fields_to_plot_mT
    tabcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, ax = plt.subplots(figsize=(6,3))
    for idx, field_mT in enumerate(fields_to_plot_mT):
        bsms_target = central_bsms + field_mT / mT_per_bsms
        kk = np.argmin(np.abs(bsms_field - bsms_target))
        bsms_true = spectra[kk]["bsms"]
        field_mT_true = (bsms_true-central_bsms) * mT_per_bsms
        ax.plot(spectra[kk]["ppm"], np.real(spectra[kk]["spec"]),
                label=f"$\\Delta B$ = {field_mT_true:.3g} mT",
                lw=2)
        ax_int.axhline(field_mT_true, color=tabcolors[idx], ls="--", alpha=.8)
        ax_map.axhline(field_mT_true, color=tabcolors[idx], ls="--")

    ax.set_xlim(max(plotRange), min(plotRange))
    ax.set_xlabel(r"$\delta$ (ppm)")
    ax.set_ylabel("NMR signal [a.u.]")
    # ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3)
    ax.grid(alpha=.3)
    ax.axvline(230, ls='--', lw=1, color="grey")
    ax.axvline(10, ls='--', lw=1, color="grey")
    fig.tight_layout()


    fig_int.savefig(os.path.join(savepath, "2dMap_figure_zoom.png"))
    fig.savefig(os.path.join(savepath, "slices.png"))
    #%%