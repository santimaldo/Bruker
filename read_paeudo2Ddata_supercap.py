# -*- coding: utf-8 -*-

import os
import numpy as np
from Datos import *

# =========================================================
# USER PARAMETERS
# =========================================================

path_local = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/"
path_bruker = "300old/2026-03-13_supercaps_YP50F_LiTFSI1M/"

savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2026-03_LiTFSI1M-aq_YP50F/"
savepath_especifico = "CA_0.75V/"

nucleo = "19F"

basepath = path_local + path_bruker
savepath = savepath_local + savepath_especifico

exp_before = 20
exp_2D = 21
exp_after = 22

ppmRange = [-49, -92]




# =========================================================
# SAVE FUNCTION (UNIVERSAL)
# =========================================================
def save_spectrum(ppm, real, imag, savepath, name):

    if not os.path.exists(savepath):
        os.makedirs(savepath)

    norm = np.max(np.abs(real))
    real_norm = real / norm
    imag_norm = imag / norm

    data = np.array([ppm, real, imag, real_norm, imag_norm]).T

    header = "ppm\t real\t imag\t real_norm\t imag_norm"
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

    file = save_spectrum(ppm, real, imag, savepath, f"{nucleo}_{label}")

    print(f"[1D] saved: {file}")


# =========================================================
# 2D EXPORT (as list of 1D spectra)
# =========================================================
def export_2D(expn, ppmRange, basepath, savepath, nucleo, label):

    datos = DatosProcesados2D(f'{basepath}{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    spectra = datos.espectro  

    for i in range(spectra.size[0]):

        ppm = spectra.ppmAxis
        real = spectra.real[i]
        imag = spectra.imag[i]

        file = save_spectrum(
            ppm, real, imag,
            savepath,
            f"{nucleo}_{label}_specn-{i}"
        )

    print(f"[2D] saved {spectra.size[0]} spectra for exp {expn}")


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