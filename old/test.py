import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = path = "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/1/"
savepath = "S:/temp/"


datos = DatosCrudos(path)




# re = datos.espectro.real
# im = datos.espectro.imag
# mod = datos.espectro.abs()

# ppmAxis = datos.espectro.ppmAxis
# ppmAxis = ppmAxis[ppmAxis>100]
# re = re[0:ppmAxis.size]

# dataexport = np.array([ppmAxis, re]).T
# #np.savetxt(savepath+'tp_'+str(pdata)+'us.dat', dataexport)

# print('graficando...')
# plt.figure()
# plt.plot(ppmAxis, re,'b')
# #plt.plot(ppmAxis, im,'r')
# #plt.plot(mod, 'g')
# plt.show()
