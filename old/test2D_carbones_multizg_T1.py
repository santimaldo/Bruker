import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *


path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY/"
savepath = "S:/temp/"

nexp = [7]
mT = [1, 10, 50, 5, 100, 25]


for j in range(len(nexp)):
    datos = DatosProcesados2D(path+str(nexp[j])+"/")
    datos.espectro.ppmSelect2D([-8, 8])
        
    ppm_x = datos.espectro.ppmGrid_Dir
    ppm_y = datos.espectro.ppmGrid_Ind
    spec = datos.espectro.real
    
    archivo_out = "mixTime_"+str(mT[j])+"ms"
    np.savetxt(savepath+archivo_out+'.dat', spec)
    np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
    np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)



##%%
#plt.figure(4563)
#plt.contourf(ppm_x, ppm_y, spec)
#plt.show()
#
#
##%%
#from mpl_toolkits.mplot3d import Axes3D
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(ppm_x, ppm_y, spec, cmap='jet')
