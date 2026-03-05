# -*- coding: utf-8 -*-
"""
Iterar y guardar espectros 1D en un rango de ppm específico
"""
import numpy as np
import matplotlib.pyplot as plt
from Datos import *

# ------------------ Configuración ------------------
expn = 62
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-06_PEO-solid-electrolyte/{expn}/"
muestra = "1H_PEO-LiTFSI"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-02_500MHz_13C-CP\PEO-LiTFSI\1H_T1/"

# Guardar archivos
save = True

# Rango de ppm que quieres seleccionar
plotRange = [20, -15]

# ------------------ Carga de datos ------------------
datos = DatosProcesadosT1(path)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spectra = datos.espectro.real  # numpy array con shape (n_spectra, n_points)
tau = datos.tau / 1000 # paso a segundos

# ------------------ Iterar y guardar ------------------
for i, spec in enumerate(spectra):
    data_to_save = np.array([ppmAxis, spec]).T  # columnas: ppm y señal
    if save:
        filename = f"{savepath}/{i}.dat"
        np.savetxt(filename, data_to_save, header="ppm\tIntensity")
        
# Guardar tau en un archivo aparte
if save:
    tau_file = f"{savepath}/tau.dat"
    np.savetxt(tau_file, tau, header="tau [s]")

# ------------------ Graficar todos los espectros ------------------
fig, ax = plt.subplots(figsize=(8,5))
cmap = plt.get_cmap("plasma")
n_spectra = spectra.shape[0]

for i, spec in enumerate(spectra):
    color = cmap(i / n_spectra)
    ax.plot(ppmAxis, spec, color=color, label=f"Spectrum {i}")

ax.set_xlim(np.max(ppmAxis), np.min(ppmAxis))  # invertir eje ppm
ax.set_xlabel("ppm")
ax.set_ylabel("Intensity")
ax.set_title(muestra)
ax.legend(fontsize=6, ncol=2)
plt.tight_layout()
plt.show()
# %%
