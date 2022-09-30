import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Exsy import *

from scipy.signal import savgol_filter


path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY/31/"


datos = DatosProcesados2D(path)

datos.espectro.ppmSelect2D([-8, 8])
# anulo todos los valores menores a 0:
datos.espectro.set_mask(0)

exsy = Exsy(datos.espectro, 3)


exsy.establecer_region((0,1),(4.5,1.9))
exsy.establecer_region((1,0),(1.9,4.5))

exsy.establecer_region((1,2),(1.9,-3.0))
exsy.establecer_region((2,1),(-3.0,1.9))

exsy.establecer_region((0,2),(4.5,-3.0))
exsy.establecer_region((2,0),(-3.0,4.5))




exsy.graficar_regiones(1)
integrales = exsy.integrar_regiones()

