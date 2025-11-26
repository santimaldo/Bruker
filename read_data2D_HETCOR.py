# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022

@author: Santi

this is for TRUE 2D data with 1D projections
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate
plt.rcParams.update({'font.size': 14})


# directorio de datos
expn = 105
expn_dir = 70
expn_ind = 82 
path  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/{expn}/"
path_dir  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/{expn_dir}/"
path_ind  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/{expn_ind}/"

# directorio de guardado
savepath= r"/"
muestra = ""

save = False
# rango de integracion
ppmRange, ppmRange_plot = [20,-20], [10, -8]
ppmRangeInd, ppmRangeInd_plot = [20,-20], [8, -4]


# === Load data ===
datos = DatosProcesados2D(path)
datos.espectro.ppmSelect2D(ppmRange, rango_ind=ppmRangeInd)
ppmAxis = datos.espectro.ppmAxis
ppmAxisInd = datos.espectro.ppmAxisInd
spec = datos.espectro.real

datos_dir = DatosProcesados(path_dir)
datos_dir.espectro.ppmSelect(ppmRange)
ppmAxis_dir = datos_dir.espectro.ppmAxis
spec_dir = datos_dir.espectro.real

datos_ind = DatosProcesados(path_ind)
datos_ind.espectro.ppmSelect(ppmRangeInd)
ppmAxis_ind = datos_ind.espectro.ppmAxis
spec_ind = datos_ind.espectro.real

# === Prepare 2D data ===
vmax = np.max(spec)
Ncontour = 10
cmap = 'grey_r'

# === Create figure with projections ===
fig = plt.figure(figsize=(7,7))
fig.suptitle(muestra)

# Define layout: 2D in center, 1D projections top/right
from matplotlib.gridspec import GridSpec
gs = GridSpec(2, 2, width_ratios=[4,1.2], height_ratios=[1.2,4], hspace=0.05, wspace=0.05)

ax_projH = fig.add_subplot(gs[0,0])   # top projection (direct, ^1H)
ax_projLi = fig.add_subplot(gs[1,1])  # right projection (indirect, ^7Li)
ax_main = fig.add_subplot(gs[1,0], sharex=ax_projH, sharey=ax_projLi)  # main 2D map

# === 2D contour plot ===
ax_main.contour(ppmAxis, ppmAxisInd, spec, Ncontour, vmax=vmax, cmap=cmap)
ax_main.set_xlim(np.max(ppmRange_plot), np.min(ppmRange_plot))
ax_main.set_ylim(np.max(ppmRangeInd_plot), np.min(ppmRangeInd_plot))
ax_main.set_xlabel(r"$^1$H Chemical Shift [ppm]")
ax_main.set_ylabel(r"$^7$Li Chemical Shift [ppm]")

# === Set ticks every 2 ppm including 0 ===
x_start = np.floor(ppmRange_plot[1]/2)*2   # nearest multiple of 2
x_end = np.ceil(ppmRange_plot[0]/2)*2
y_start = np.floor(ppmRangeInd_plot[1]/2)*2
y_end = np.ceil(ppmRangeInd_plot[0]/2)*2
ax_main.set_xticks(np.arange(x_start, x_end + 1, 2))
ax_main.set_yticks(np.arange(y_start, y_end + 1, 2))
ax_main.grid(True,  ls='--', alpha=0.8)


# === 1D projections ===
# scale for better visual alignment
scaleH = 0.3 * np.max(ppmAxisInd)  # vertical scaling
scaleLi = 0.3 * np.max(ppmAxis)    # horizontal scaling

ax_projH.plot(ppmAxis_dir, spec_dir/np.max(spec_dir)*np.max(ppmAxisInd)*0.2, color='k')
ax_projH.set_xlim(np.max(ppmRange_plot), np.min(ppmRange_plot))
ax_projH.axis('off')

spec_ind = np.sum(spec, axis=1)  ## sum along direct dimension
ppmAxis_ind = ppmAxisInd
ax_projLi.plot(spec_ind/np.max(spec_ind)*np.max(ppmAxis)*0.2 + np.min(ppmAxisInd), ppmAxis_ind, color='k')
ax_projLi.set_ylim(np.max(ppmRangeInd_plot), np.min(ppmRangeInd_plot))
ax_projLi.axis('off')



#%% === Slices along indirect dimension ===
# Horizontal slices: fixed ^7Li ppm (indirect dimension)
slice_ppm_values_ind = [1, 3.1]   # <-- cámbialos según tus regiones de interés

fig, ax = plt.subplots(1, 1, figsize=(8,6))
colors = ['C0', 'C1']

for i, ppm_val in enumerate(slice_ppm_values_ind):
    idx = np.argmin(np.abs(ppmAxisInd - ppm_val))
    slice_spec = spec[idx, :]
    slice_spec = slice_spec / np.max(np.abs(slice_spec))
    ax.plot(ppmAxis, slice_spec, color=colors[i], lw=3,
            label=f"{ppmAxisInd[idx]:.1f} ppm "+r"($^7$Li)")
    
ax.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax.set_xlabel(r"$^1$H Chemical Shift [ppm]")
ax.set_ylabel("Normalized Intensity (offset for clarity)")
ax.legend()
ax.set_title(r"Slices along indirect dimension (fixed $^7$Li ppm)")
plt.tight_layout()


plt.show()


#%% === Slices along direct dimension ===
# Vertical slices: fixed ^1H ppm (direct dimension)
slice_ppm_values_dir = [-1, 4]   # <-- cámbialos según tus regiones de interés

fig, ax = plt.subplots(1, 1, figsize=(8,6))
colors = ['C2', 'C3']

for i, ppm_val in enumerate(slice_ppm_values_dir):
    idx = np.argmin(np.abs(ppmAxis - ppm_val))
    # extract vertical slice (i.e. along indirect dimension)
    slice_spec = spec[:, idx]
    slice_spec = slice_spec / np.max(np.abs(slice_spec))
    ax.plot(ppmAxisInd, slice_spec, color=colors[i],
            label=f"{ppmAxis[idx]:.1f} ppm ($^1$H)")
    
ax.set_xlim(np.max(ppmAxisInd), np.min(ppmAxisInd))
ax.set_xlabel(r"$^7$Li Chemical Shift [ppm]")
ax.set_ylabel("Normalized Intensity (offset for clarity)")
ax.legend()
ax.set_title(r"Slices along direct dimension (fixed $^1$H ppm)")

plt.tight_layout()
plt.show()





# === Save data if needed ===
if save:
    filename = f'{savepath}/{muestra}.png'
    fig.savefig(filename, dpi=300)
    np.savetxt(f"{savepath}/{muestra}_data2D.dat", spec)
    np.savetxt(f"{savepath}/{muestra}_ppmAxis.dat", ppmAxis)
    np.savetxt(f"{savepath}/{muestra}_ppmAxisInd.dat", ppmAxisInd)

plt.show()

# %%
