# -*- coding: utf-8 -*-

import os
import numpy as np
from Datos import *

# =========================================================
# USER PARAMETERS
# =========================================================

# ########################################## SUPERCAP 2026/03
# path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/"
# path_bruker = "300old/2026-03-13_supercaps_YP50F_LiTFSI1M/"

# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2026-03_LiTFSI1M-aq_YP50F/"
# # savepath_especifico = "CA_-0.75V/"
# # exp_before = 20
# # exp_2D = 21
# # exp_after = 22
# ###--------------------------------
# # savepath_especifico = "CA_0.75V/"
# # exp_before = 30
# # exp_2D = 31
# # exp_after = 32
# ###--------------------------------
# savepath_especifico = "AC_10Hz/"
# exp_before = 40
# exp_2D = 41
# exp_after = 42

# nucleo = "19F"
# basepath = path_local + path_bruker
# savepath = savepath_local + savepath_especifico
# # - -  - - - - - - - - - - - - - - - - - - - - - -
# ppmRange = [-50, -92]

########################################## SUPERCAP 2026/03
path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/"
path_bruker = "300old/2026-04-07_supercaps_YP50F_LiTFSI1M/"

savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2026-04_LiTFSI1M-aq_YP50F/"
###--------------------------------
savepath_especifico = "CA_-0.25V_expn-61/"
exp_before = 60
# ###--------------------------------
# savepath_especifico = "CA_0.25V/"
# exp_before = 70
# ###--------------------------------
# savepath_especifico = "CA_0.75V/"
# exp_before = 80
# ###--------------------------------
# savepath_especifico = "CA_-0.50V/"
# exp_before = 90
# ###--------------------------------
# savepath_especifico = "CA_0.50V/"
# exp_before = 100
# ###--------------------------------
# savepath_especifico = "CA_-0.25V_expn-111/"
# exp_before = 110
# ###--------------------------------
# savepath_especifico = "CA_1.00V/"
# exp_before = 120
# ###--------------------------------
# savepath_especifico = "CA_-1.00V/"
# exp_before = 130
###--------------------------------
savepath_especifico = "CA_-0.75V/"
exp_before = 140

exp_2D = exp_before + 1
exp_after = exp_before + 2
nucleo = "19F"
basepath = path_local + path_bruker
savepath = savepath_local + savepath_especifico
# - -  - - - - - - - - - - - - - - - - - - - - - -
ppmRange = [-50, -92]




# =========================================================
# SAVE FUNCTION (UNIVERSAL)
# =========================================================
def save_spectrum(ppm, real, imag, savepath, name):

    if not os.path.exists(savepath):
        os.makedirs(savepath)    

    data = np.array([ppm, real, imag]).T

    header = "ppm\t real\t imag"
    filename = os.path.join(savepath, f"{name}.dat")

    np.savetxt(filename, data, header=header)

    return filename


# =========================================================
# 1D EXPORT (before / after)
# =========================================================
def export_1D(expn, ppmRange, label, basepath, savepath, nucleo):

    datos = DatosProcesados(f'{basepath}{expn}/')
    datos.espectro.ppmSelect(ppmRange)
    ppm = datos.espectro.ppmAxis
    real = datos.espectro.real
    imag = datos.espectro.imag

    fig, ax = plt.subplots()
    ax.plot(ppm, real, label=f"expn: {expn}")
    ax.set_xlabel("ppm")
    ax.set_ylabel("Intensity")
    ax.legend()

    file = save_spectrum(ppm, real, imag, savepath, f"{nucleo}_{label}")

    print(f"[1D] saved: {file}")


# =========================================================
# 2D EXPORT (as list of 1D spectra)
# =========================================================
def export_2D(expn, ppmRange, basepath, savepath, nucleo, label):

    datos = DatosProcesados2D(f'{basepath}{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    spectra = datos.espectro  
    Number_valid_data = int(np.where(spectra.real[:, 0]==0)[0][0]) - 1
    spec_and_time = []
    for i in range(Number_valid_data+1):

        ppm = spectra.ppmAxis
        real = spectra.real[i]
        imag = spectra.imag[i]
        if real.max() - real.min() < 1e-6:
            continue  # Skip spectra with 0 signal
        file = save_spectrum(
            ppm, real, imag,
            savepath,   
            f"{nucleo}_{label}_specn-{i}"
        )
    
    Info = np.array([Number_valid_data, datos.acqus.D1])
    np.savetxt(os.path.join(savepath, f"Info"), Info, header="Number of valid spectra\tD1(s)")

    print(f"[2D] saved {Number_valid_data + 1} spectra for exp {expn}")


# =========================================================
# RUN PIPELINE
# =========================================================

# -------------------
# 1D BEFORE
# -------------------
export_1D(
    expn=exp_before,
    ppmRange=ppmRange,
    label="1d-before",
    basepath=basepath,
    savepath=savepath,
    nucleo=nucleo
)

# -------------------
# 2D
# -------------------
export_2D(
    expn=exp_2D,
    ppmRange=ppmRange,
    basepath=basepath,
    savepath=savepath,
    nucleo=nucleo,
    label="2D"
)

# -------------------
# 1D AFTER
# -------------------
export_1D(
    expn=exp_after,
    ppmRange=ppmRange,
    label="1d-after",
    basepath=basepath,
    savepath=savepath,
    nucleo=nucleo
)