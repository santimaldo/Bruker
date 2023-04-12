import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

from scipy.signal import savgol_filter


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]



path  = "S:/Doctorado/Carbones/300MHz/2021-12-23_Carbones_MAS/18/"
savepath = "S:/tmp/"
muestra = "SAC-1H-1ms"


datos = DatosProcesados2D(path)

datos.espectro.ppmSelect2D([-15,6])


ppm_x = datos.espectro.ppmGridDir
ppm_y = datos.espectro.ppmGridInd


ppmDir = datos.espectro.ppmAxis
ppmInd = datos.espectro.ppmAxisInd



spec = datos.espectro.real

mask = 0 
spec[spec<mask] = 0
#%%

plt.figure(4563)
plt.contour(ppm_x, ppm_y, spec, 100,  cmap='jet', vmax=5000000)
ax = plt.gca()
ax.invert_yaxis()
ax.invert_xaxis()
plt.show()




np.savetxt(f"{savepath}{muestra}_data2D.dat", spec)
np.savetxt(f"{savepath}{muestra}_ppmDir.dat", ppm_x)
np.savetxt(f"{savepath}{muestra}_ppmInd.dat", ppm_y)