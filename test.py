import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/1/"


datos = DatosProcesados1D(path)
datos.set_espectro()



print(datos.acqus.NS)


re = datos.espectro.real
im = datos.espectro.imag
mod = datos.espectro.abs()


print('graficando...')
plt.figure()
plt.plot(re,'b')
plt.plot(im,'r')
plt.plot(mod, 'g')
plt.show()
