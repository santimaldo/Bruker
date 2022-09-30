import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  =  "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/"
savepath = "S:/Doctorado/Carbones/analisis/2020-01_CM7_EXSY/muestra2019+agua/"

nexp = [1019]

for j in range(len(nexp)):
    datos = DatosProcesados2D(path+str(nexp[j])+"/")
    # datos.espectro.ppmSelect2D([-8, 8])
  
    spec = datos.espectro.real
    
    ppm_Axis = datos.espectro.ppmAxis
    
##%%
tp = np.arange(10,321,10)
plt.figure(4563)
plt.contourf(ppm_Axis, tp , spec)
plt.show()
#
#
# ##%%
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure(4568)
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(tp, ppm_Axis , spec, cmap='jet')
