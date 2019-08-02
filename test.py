import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

# en $olvia
#path = "../2019-07-31_Esferas_CantidadDeHumedad/3/"
# en $lily
path = "S:/Doctorado/Dendritas/116MHz/2019-06-13_Dendritas/1/"



datos = DatosProcesados(path)

print(datos.acqus.NS)

re = datos.espectro.real
im = datos.espectro.imag
mod = datos.espectro.abs()

print('graficando...')
plt.figure()
plt.plot(re, 'b')
plt.plot(im, 'r')
plt.plot(mod, 'g')
plt.show()
