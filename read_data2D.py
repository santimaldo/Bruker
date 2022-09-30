# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate




# directorio de datos
expn = 122
path  =f"S:/CarbonesSofi/300MHz/2022-05-12_Carbones_CMK3ACT/{expn}/" # compu Ofi
# directorio de guradado
savepath= "S:/CarbonesSofi/Analisis/2022-09_EXSY/data2D/"
muestra = "1H_EXSY_CMK3ACT_mT_001ms"


save = False
#rango de integracion
ppmRange = [50,-50]

datos = DatosProcesados2D(path)
datos.espectro.ppmSelect2D(ppmRange)
ppmAxis = datos.espectro.ppmAxis
ppmAxisInd = datos.espectro.ppmAxisInd
spec = datos.espectro.real


# spec[spec<  1/100*np.max(spec)] = 0
vmax = 0.8*np.max(spec)
Ncontour = 80
cmap = 'inferno'

fig, ax = plt.subplots(1,1)
fig.suptitle(muestra)
ax.contour(ppmAxis, ppmAxisInd, spec,Ncontour, vmax=vmax, cmap=cmap)
ax.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
ax.set_ylim(np.max(ppmAxisInd), np.min(ppmAxisInd))
ax.set_xlabel(r"$^1$H Chemical Shift [ppm]")
ax.set_ylabel(r"$^1$H Chemical Shift [ppm]")


# guardo data:
if save:
  filename = f'{savepath}/{muestra}.png'
  fig.savefig(filename)   # save the figure to file


  np.savetxt(f"{savepath}/{muestra}_data2D.dat", spec)
  np.savetxt(f"{savepath}/{muestra}_ppmAxis.dat", ppmAxis)
  np.savetxt(f"{savepath}/{muestra}_ppmAxisInd.dat", ppmAxisInd)







