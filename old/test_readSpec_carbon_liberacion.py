import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *



# directorio de guradado
path1 = "S:/Doctorado/Carbones/analisis/2019-08_Liberacion/CM7/Remedicion/spec/"
path1 = "S:/temp/"

# directorio de datos
path0 = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"


tiempos = [0,0.5,1,2,4,8,13,40,61,80,170]
expnum = [i*4+1 for i in range(len(tiempos))] # litio
expnum = [i*4+1+2 for i in range(len(tiempos))] # agua
integrales = []
integrales_fit = []

for expn in expnum:    
    directorio = path0+str(expn)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    
    spec = datos.espectro.real / NS    
    ppmAxis = datos.espectro.ppmAxis
#    ppmAxis = ppmAxis[ppmAxis>100]
#    re = re[0:ppmAxis.size]        
    vf = VoigtFit(ppmAxis,spec, Npicos=2)
    ajuste, componentes = vf.componentes(ppmAxis)
    
    integrales.append(integrate.simps(spec, x=ppmAxis))
    integrales_fit.append(integrate.trapz(ajuste))
    
#    archivo_out = 't_170h'
#    dataexport = np.array([ppmAxis, re]).T
    #np.savetxt(savepath+'tp_'+str(pdata)+'us.dat', dataexport)
    
    print('graficando...')
    plt.figure(expn)
    plt.plot(ppmAxis, spec)
    plt.plot(ppmAxis, ajuste, 'k')
    for comp in componentes:
        plt.plot(ppmAxis, comp, '--')

#%%   
        
# # espectro sobre el cual normalizar        
# directorio = path0+str(1)+"/"
# datos = DatosProcesados(directorio)    
# NS = datos.acqus.NS  

# spec = datos.espectro.real / NS    
# ppmAxis = datos.espectro.ppmAxis
# #    ppmAxis = ppmAxis[ppmAxis>100]
# #    re = re[0:ppmAxis.size]        
# vf = VoigtFit(ppmAxis,spec, Npicos=2)
# ajuste, componentes = vf.componentes(ppmAxis)

# FactorNorm = integrate.simps(spec, x=ppmAxis)
# FactorNorm_fit = integrate.trapz(ajuste)

        
# plt.figure(789)
# integrales = integrales/FactorNor
# integrales_fit = integrales_fit/FactorNorm_fit
# plt.plot(tiempos, integrales,'o-', label='Integral de los datos')
# plt.plot(tiempos, integrales_fit,'o-', label='Integral del ajuste')
# plt.show()

#%%
tiempos = np.array(tiempos)
integrales = np.array(integrales)
media = np.mean(integrales)
error = 2*np.std(integrales)/np.sqrt(integrales.size)

integrales=integrales/media
error = error/media
media = 1


plt.figure(78756)

plt.plot(tiempos, integrales,'ro-', label='Integral de los datos')
plt.plot(tiempos, tiempos*0+media,'k-', label='media')
plt.plot(tiempos, tiempos*0+media+error,'k--', label='media')
plt.plot(tiempos, tiempos*0+media-error,'k--', label='media')
plt.xlabel("Tiempo [h]")

plt.show()
