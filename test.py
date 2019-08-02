import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = "../2019-07-31_Esferas_CantidadDeHumedad/3/"


datos = DatosProcesados(path)



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
