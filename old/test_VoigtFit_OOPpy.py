# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:37:40 2019

@author: santi
"""
import numpy as np
import matplotlib.pyplot as plt
from VoigtFit import *
from Datos import *


# directorio de datos
path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"
path  = "S:/Doctorado/LiMetal/116MHz/2019-09-11_Dendritas_SiO2/"
expn = 2
directorio = path+str(expn)+"/"
datos = DatosProcesados(directorio)
NS = datos.acqus.NS
spec = datos.espectro.real / NS
ppmAxis = datos.espectro.ppmAxis
ppmAxis = ppmAxis[ppmAxis>100]
spec = spec[0:ppmAxis.size]

# vf = VoigtFit(ppmAxis,spec, Npicos=3, center=[244,252,260], fijar=['m3_center'])
vf = VoigtFit(ppmAxis,spec, Npicos=3, center=[261,244])
ajuste, componentes = vf.componentes(ppmAxis)


plt.figure(1231)
plt.plot(ppmAxis,spec)
plt.plot(ppmAxis, ajuste, 'k')
for comp in componentes:
    plt.plot(ppmAxis, comp, '--')
# plt.xlim([350,150])


# params = vf.params
# vf = VoigtFit(ppmAxis,spec, Npicos=2, params=params)
# ajuste, componentes = vf.componentes(ppmAxis)


# plt.figure(1232)
# plt.plot(ppmAxis,spec)
# plt.plot(ppmAxis, ajuste, 'k')
# for comp in componentes:
#     plt.plot(ppmAxis, comp, '--')
# plt.xlim([350,150])


plt.show()