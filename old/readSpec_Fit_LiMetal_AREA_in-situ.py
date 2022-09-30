import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *



# directorio de guradado

savepath = "S:/Doctorado/LiMetal/Analisis/2021-08_in-situ/celda3/"

# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2021-08-18_litio_in-situ/"



# tiempos = [0,1,4,7,13,21,31,45,60,75,90,120,150,180,210] # celda 1
# tiempos = [0,5,10,20,30,45,60,90] # celda 2
tiempos = [0,5,10,20,30,45,60,65,75,90,120] # celda 3

expnum = np.arange(1,len(tiempos)+1) +1
integrales = []
integrales_fit = []
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)
    
    NS = datos.acqus.NS
    RG = datos.acqus.RG
    
    spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    spec = spec[ppmAxis>100]
    ppmAxis = ppmAxis[ppmAxis>100]
    spec = spec[ppmAxis<400]
    ppmAxis = ppmAxis[ppmAxis<400]
    
  
    # vf = VoigtFit(ppmAxis,spec, Npicos=3)
    # ajuste, componentes = vf.componentes(ppmAxis)
    
    integrales.append(integrate.trapz(spec))
    # integrales_fit.append(integrate.trapz(ajuste))
    
    archivo_out = 'espectro'+str(expn)+'.dat'
    dataexport = np.array([ppmAxis, spec]).T
    np.savetxt(savepath+archivo_out, dataexport)
    
    print('graficando...')
    plt.figure(123568)
    plt.plot(ppmAxis, spec)
    # plt.plot(ppmAxis, ajuste, 'k')
    # for comp in componentes:
    #     plt.plot(ppmAxis, comp, '--')

#%%   
plt.figure(789)
integrales = integrales/integrales[0]
# integrales_fit = integrales_fit/max(integrales_fit)
plt.plot(tiempos, integrales,'o', label='Integral de los datos')
# plt.plot(tiempos, integrales_fit,'o-', label='Integral del ajuste')
plt.show()


archivo_out = 'integrales.dat'
dataexport = np.array([tiempos, integrales]).T
np.savetxt(savepath+archivo_out, dataexport)