# -*- coding: utf-8 -*-
"""
Created on Feb 20/02/2025


@author: Santi
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
from Datos import *
plt.rcParams['font.size'] = 14

expn_before = 70    # Li
expn_pseudo2d = 71  # Li


# Data paths
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/300old/2025-02-07_insitu-sync-start/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-02_LiTFSI1M-aq_CA-cycles/"
muestra = "7Li_"
save = False
plotRange = [4, -8]

# Define colormaps and experiments
colormap = "Reds_r"

# Load 1D experiment data
datos = DatosProcesados(f'{path}/{expn_before}/')
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
re1d = datos.espectro.real

# Load 2D experiment data
datos = DatosProcesados2D(f'{path}/{expn_pseudo2d}/')
datos.espectro.ppmSelect(plotRange)
vdlist = datos.get_vdlist() / 1000  # Convert to seconds
datas = [datos.espectro.real[:vdlist.size], datos.espectro.real[vdlist.size:]]

titles = [r"0 V $\rightarrow$ 1 V", r"1 V $\rightarrow$ 0 V"]
# Create figure with two subplots (stacked vertically)
fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=False)
for ax in axs:
    ax.axhline(y=0, color='grey', linestyle='-')
    ax.set_xlabel(r"$^7$Li $\Delta\delta$ [ppm]")
for idx, (data, ax) in enumerate(zip(datas, axs)):

    num_curves = vdlist.size
    cmap = cm.get_cmap(colormap, num_curves)  # Discrete colormap
    bounds = np.arange(num_curves + 1) - 0.5  # Boundaries
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)

    colors = [cmap(i) for i in np.linspace(0.05, 0.6, num_curves)]

    ax.plot(ppmAxis, re1d, '--', color="black", label="Before pseudo-2D")
    # Plot each spectrum in its corresponding subplot
    for i, (re, color) in enumerate(zip(data, colors)):
        ax.plot(ppmAxis, re, color=color)

    ax.set_xlim(plotRange)
    ax.text(0.05, 0.9, titles[idx], transform=ax.transAxes, fontsize=16, verticalalignment='top')
    ax.legend()
    # Add colorbar to each subplot
    cbar = plt.colorbar(sm, ax=ax, ticks=np.arange(num_curves))
    cbar.set_label("Variable delay [s]")
    cbar.set_ticks(np.arange(num_curves))
    cbar.set_ticklabels([f"{vdlist[i]}" for i in range(num_curves)], fontsize=10)

# Set common x-axis label
axs[-1].set_xlabel(r"$^7$Li $\Delta\delta$ [ppm]")

plt.tight_layout()  # Adjust spacing between subplots
plt.show()





# # guardo data:
# if save:
#     filename = f'{savepath}/{muestra}_T1.png'
#     fig.savefig(filename)   # save the figure to file

#     Signals = np.array(Signals).T
#     vdlist = vdlist.reshape(vdlist.size, 1)
#     T1data = np.hstack((vdlist, Signals))
#     header = "vdlist [s]\t"
#     for ppmRange in ppmRanges:
#         header += f"{ppmRange} ppm\t"
#     np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

#     data = np.array([ppmAxis, re, im]).T
#     np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)
