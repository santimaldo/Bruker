import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1022/"
datos = DatosProcesados(path)


ppmAxis = datos.espectro.ppmAxis
nptsT1 = datos.espectro.size[0]

plt.figure(4563)
for j in range(nptsT1):    
    spec = datos.espectro.real[j]
    plt.plot(ppmAxis, spec)    

plt.figure(4564)
for j in range(nptsT1):    
    spec = datos.espectro.imag[j]
    plt.plot(ppmAxis, spec)    

plt.show()
