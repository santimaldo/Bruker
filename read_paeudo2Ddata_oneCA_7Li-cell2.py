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

# 0 V to 1 V
expn_pos_before = 4    # Li
expn_pos_pseudo2d = 6 # Li
# 0 V to -1 V
expn_neg_before = 8  # Li
expn_neg_pseudo2d = 9  # Li

# Data paths
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata/300old/2025-03-10_insitu-LiTFSIaq-supercap/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY/oneCA/"
muestra = "7Li_"
save = False
plotRange = [4, -8]

# Define colormaps and experiments
colormaps = ["Reds_r", "Blues_r"]
expn_before_list = [expn_pos_before, expn_neg_before]
expn_pseudo2d_list = [expn_pos_pseudo2d, expn_neg_pseudo2d]
titles = [r"0 V $\rightarrow$ 1 V", r"0 V $\rightarrow$ -1 V"]

# Create figure with two subplots (stacked vertically)
fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=False)
for ax in axs:
    ax.axhline(y=0, color='grey', linestyle='-')
    ax.set_xlabel(r"$^7$Li $\Delta\delta$ [ppm]")
for idx, (colormap, expn_before, expn_pseudo2d, ax) in enumerate(zip(colormaps, expn_before_list, expn_pseudo2d_list, axs)):
    # Load 1D experiment data
    datos = DatosProcesados(f'{path}/{expn_before}/')
    datos.espectro.ppmSelect(plotRange)
    ppmAxis1d = datos.espectro.ppmAxis
    re1d = datos.espectro.real
    
    # Load 2D experiment data
    datos = DatosProcesados2D(f'{path}/{expn_pseudo2d}/')
    datos.espectro.ppmSelect(plotRange)
    spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    vdlist = datos.acqus.D1 * np.arange(datos.acqu2s.TD)  # Convert to seconds

    # Define colormap and colors for each curve
    num_curves = vdlist.size
    cmap = cm.get_cmap(colormap)  # Get full colormap
    color_range = np.linspace(0.05, 0.7, num_curves)  # Define the specific range of the colormap to use
    colors = [cmap(i) for i in color_range]  # Extract the exact colors used for plotting
    # Create a discrete colormap matching the plotted colors
    discrete_cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(num_curves + 1) - 0.5
    norm = mcolors.BoundaryNorm(bounds, num_curves)
    sm = cm.ScalarMappable(cmap=discrete_cmap, norm=norm)

    colors = [cmap(i) for i in np.linspace(0.05, 0.6, num_curves)]

    # Plot each spectrum in its corresponding subplot
    for i, (re, color) in enumerate(zip(datos.espectro.real, colors)):
        ax.plot(ppmAxis, re, color=color)
    ax.plot(ppmAxis1d, re1d, '--', color="black", label="Before pseudo-2D")

    ax.set_xlim(plotRange)
    ax.text(0.05, 0.9, titles[idx], transform=ax.transAxes, fontsize=16, verticalalignment='top')
    ax.legend()
    # Add colorbar to each subplot
    cbar = plt.colorbar(sm, ax=ax, ticks=np.arange(num_curves))
    cbar.set_label("Variable delay [s]")
    cbar.set_ticks(np.arange(num_curves))
    cbar.set_ticklabels([f"{vdlist[i]}" for i in range(num_curves)], fontsize=10)



plt.tight_layout()  # Adjust spacing between subplots
plt.show()

stop



# guardo data:
if save:
    filename = f'{savepath}/{muestra}_T1.png'
    fig.savefig(filename)   # save the figure to file

    Signals = np.array(Signals).T
    vdlist = vdlist.reshape(vdlist.size, 1)
    T1data = np.hstack((vdlist, Signals))
    header = "vdlist [s]\t"
    for ppmRange in ppmRanges:
        header += f"{ppmRange} ppm\t"
    np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data, header=header)

    data = np.array([ppmAxis, re, im]).T
    np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)
