import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *



# directorio de guradado
path1 = "S:/Doctorado/Carbones/analisis/2019-08_Liberacion/CM7/Remedicion/spec/"

# directorio de datos
path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"

tiempos = [0,1,4,13,27]
expnum = [i*4+1 for i in range(len(tiempos))]
integrales = []
integrales_fit = []
for expn in expnum:
    
    
        
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)
    
    NS = datos.acqus.NS
    
    
    spec = datos.espectro.real / NS
    
    ppmAxis = datos.espectro.ppmAxis
#    ppmAxis = ppmAxis[ppmAxis>100]
#    re = re[0:ppmAxis.size]
        
    vf = VoigtFit(ppmAxis,spec, Npicos=2)
    ajuste, componentes = vf.componentes(ppmAxis)
    
    integrales.append(integrate.trapz(spec))
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
plt.figure(789)
integrales = integrales/max(integrales)
integrales_fit = integrales_fit/max(integrales_fit)
plt.plot(tiempos, integrales,'o-', label='Integral de los datos')
plt.plot(tiempos, integrales_fit,'o-', label='Integral del ajuste')
plt.show()