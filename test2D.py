import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/34/"

datos = DatosProcesados(path)


ppmAxis = datos.espectro.ppmAxis
nptsT1 = datos.espectro.size[0]

plt.figure(4563)
for j in range(nptsT1):    
    spec = datos.espectro.real[j]
    plt.plot(ppmAxis, spec)    

plt.show()
