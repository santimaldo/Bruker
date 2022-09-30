import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *


path  ="S:/Doctorado/Carbones/300MHz/2021-12-27_Carbones_MAS/"
savepath = "S:/temp/"


mT =   [200]
nexp = [21 ]


I32 = []
I23 = []
plt.figure(1)
for j in range(len(nexp)):
    print(path+str(nexp[j])+"/", "  mT=", str(mT[j]))
    datos = DatosProcesados2D(path+str(nexp[j])+"/")
    datos.espectro.ppmSelect2D([-20, 20])
        
    ppm_x = datos.espectro.ppmGrid_Dir
    ppm_y = datos.espectro.ppmGrid_Ind
    spec = datos.espectro.real
    
#    archivo_out = "mixTime_"+str(mT[j])+"ms"
#    np.savetxt(savepath+archivo_out+'.dat', spec)
#    np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
#    np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)


    plt.subplot(2,3,j+1)
    plt.title('regiones a integrar')
    plt.contour(ppm_x, ppm_y, spec, 20,  cmap='jet', vmax=5000000)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.invert_xaxis()
    #slices x=-3.0 y=1.9
    slice_x = ppm_x[yi_32:yf_32,xi_32:xf_32]
    slice_y = ppm_y[yi_32:yf_32,xi_32:xf_32]
    slice_spec = spec[yi_32:yf_32,xi_32:xf_32]
    plt.contourf(slice_x, slice_y, slice_spec,  cmap='jet')
    slice_ppmDir = ppmDir[xi_32:xf_32]
    slice_ppmInd = ppmInd[yi_32:yf_32]
    I32.append(integrar(slice_ppmDir, slice_ppmInd, slice_spec))
    
    #slices x=1.9, y=-3.0 
    slice_x = ppm_x[yi_23:yf_23,xi_23:xf_23]
    slice_y = ppm_y[yi_23:yf_23,xi_23:xf_23]
    slice_spec = spec[yi_23:yf_23,xi_23:xf_23]
    plt.contourf(slice_x, slice_y, slice_spec,  cmap='jet')
    slice_ppmDir = ppmDir[xi_23:xf_23]
    slice_ppmInd = ppmInd[yi_23:yf_23]
    I23.append(integrar(slice_ppmDir, slice_ppmInd, slice_spec))


plt.figure(2)
plt.plot(mT,I32,'o')
plt.plot(mT,I23,'o')

plt.show()


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
