# -*- coding: utf-8 -*-

import os
import numpy as np
from Datos import *

# =========================================================
# USER PARAMETERS
# =========================================================

########################################## SUPERCAP 2026/03
path_local = r"C:\Users\Santi\Documents\NMRdata\400dnp/"

path_bruker = "20260716_Cambridge_Grey_LiMetal_KlystronTest/"

# savepath_local = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2026-04_LiTFSI1M-aq_YP50F/"
###--------------------------------
# savepath_especifico = "CA_-0.25V_expn-61/"


expns_2D = [30, 31]
expn_offOE = 24
nucleo = "7Li"
basepath = path_local + path_bruker
# savepath = savepath_local + savepath_especifico
# - -  - - - - - - - - - - - - - - - - - - - - - -
ppmRange = [380, 80]



integrals = np.array([])
times = np.array([])
spectra_to_extract=[[0, 133, 136], [313]]
spectra_extracted = []
spectra_extracted_integral = []
for ii, expn in enumerate(expns_2D):
    datos = DatosProcesados2D(f'{basepath}{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    spectra = datos.espectro     
    spec_and_time = []
    
    integral = datos.Integrar2D()
    integral = integral[integral!=0]
    integrals = np.concatenate([integrals, integral])
    start_time = datos.acqus.dic['DATE_START']
    end_time = datos.acqus.dic['DATE']
    time = np.linspace(start_time, end_time, integral.size)
    times = np.concatenate([times, time])

    spectra_extracted_ii = []
    spectra_extracted_integral_ii = []
    for jj in spectra_to_extract[ii]:
        ppm = spectra.ppmAxis
        real = spectra.real[jj]
        data_to_extraxt = np.array([ppm, real])
        spectra_extracted_ii.append(data_to_extraxt)
        spectra_extracted_integral_ii.append(integral[jj])
    spectra_extracted.append(spectra_extracted_ii)
    spectra_extracted_integral.append(spectra_extracted_integral_ii)
times -= times[0]
times_min = times/60

datos = DatosProcesados(f'{basepath}{expn_offOE}/')
datos.espectro.ppmSelect(ppmRange)
ppm_offOE = datos.espectro.ppmAxis
spec_offOE = datos.espectro.real
integral_offOE = datos.Integrar()

print("======== Enhancements ===========")
for exp in spectra_extracted_integral:
    for integral_ON in exp:
        print(integral_ON/integral_offOE)
print("=================================")
#%%
fig, ax = plt.subplots(figsize=(8,4))
ax.plot(times_min, integrals/integrals[0], 'ko-')
ax.set_xlabel("Time [min]")
ax.set_ylabel("Integral (norm. to time zero)")
ax.set_xlim(225,240)
#%%
figsize = (4,3)
fig, ax = plt.subplots(figsize=figsize)
ax.axhline(0, ls='--', lw="0.5", color='k')
ax.plot(ppm_offOE, spec_offOE, "grey")
ppm, spec = spectra_extracted[0][0]
ax.plot(ppm, spec, "tab:blue")
ppm, spec = spectra_extracted[1][0]
ax.plot(ppm, spec, color="tab:orange")
ax.set_xlim([max(ppmRange), min(ppmRange)])
ax.set_xlabel("ppm")

fig, ax = plt.subplots(figsize=figsize)
ax.axhline(0, ls='--', lw="0.5", color='k')
ax.plot(ppm_offOE, spec_offOE, "grey")
ppm, spec = spectra_extracted[0][1]
ax.plot(ppm, spec, color="tab:green")
ppm, spec = spectra_extracted[0][2]
ax.plot(ppm, spec, color="tab:red")
ax.set_xlim([max(ppmRange), min(ppmRange)])
ax.set_xlabel("ppm")


fig, ax = plt.subplots(figsize=figsize)
ax.axhline(0, ls='--', lw="0.5", color='k')
ax.plot(ppm_offOE, spec_offOE, "grey")
ppm, spec = spectra_extracted[0][0]
ax.plot(ppm, spec, "tab:blue")
ppm, spec = spectra_extracted[1][0]
ax.plot(ppm, spec, color="tab:orange")
ax.set_xlim([max(ppmRange), min(ppmRange)])
ax.set_xlabel("ppm")

ppm, spec = spectra_extracted[0][1]
ax.plot(ppm, spec, '--', color="tab:green")
ppm, spec = spectra_extracted[0][2]
ax.plot(ppm, spec, '--', color="tab:red")
ax.set_xlim([280, 200])
ax.set_xlabel("ppm")

# %%
