import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *




# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2021-11-08_cabezal_litio_senal_vs_z/"


# celda S1
alturas = [0 , 0.75, 1.5, 3 , 5 ] 
expnum  = [21, 26  , 24 , 25, 22 ]
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_B1_vs_z/Swagelok1/"

# # celda S2
# alturas = [0, 0.75, 1.5, 3, 4.5 ] 
# expnum  = [5, 6   , 4  , 2, 3 ]
# # directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_B1_vs_z/Swagelok2/"

# Aluminas
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2021-11-23_Aluminas/"
expnum  = [4, 20, 31, 34, 36]
filenames = ["LiMetal", "LiMetal_Al2O3", "Blanco_DEP", "Blanco_DIS", "Al2O3_DEP" ]
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_ALuminas/"


integrales = []
integrales_err = []
jj=0
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)
    
    # NS = datos.acqus.NS
    # RG = datos.acqus.RG
    
        
    spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    ppm_med = ppmAxis[spec==np.max(spec)][0]
    
    i = 0
    integrales_i = []
    for ppm_window in np.arange(80,201,10):
      ppm_dif = ppm_window/2
      condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
      spec_i = spec[condicion]        
      ppmAxis_i = ppmAxis[condicion]    
      
      integrales_i.append(integrate.simps(spec_i, x=ppmAxis_i))
   
    integrales.append(np.mean(integrales_i))
    integrales_err.append(np.std(integrales_i)*2)
    archivo_out = filenames[jj]+'.dat'
    dataexport = np.array([ppmAxis, spec]).T
    np.savetxt(savepath+archivo_out, dataexport)
    # print('graficando...')
    plt.figure(123568)
    plt.plot(ppmAxis, spec)
    
    jj+=1

integrales = np.array(integrales)    
integrales_err = np.array(integrales_err)
#%%   
plt.figure(789)
# integrales = integrales/integrales[0]
# integrales_fit = integrales_fit/max(integrales_fit)
# plt.errorbar(alturas, integrales, yerr=integrales_err, fmt='o')
# plt.plot(alturas, integrales,'ko', label='Integral de los datos')
# plt.plot(tiempos, integrales_fit,'o-', label='Integral del ajuste')
plt.errorbar(np.arange(integrales.size),integrales, yerr=integrales_err, fmt='o')
plt.show()


integrales_norm = integrales/integrales[0]
integrales_norm_err = integrales_norm*np.sqrt((integrales_err/integrales)**2+(integrales_err[0]/integrales[0])**2)
plt.figure(7810)
plt.errorbar(alturas, integrales_norm, yerr=integrales_norm_err, fmt='o')



archivo_out = 'integrales.dat'
dataexport = np.array([alturas, integrales, integrales_err,integrales_norm, integrales_norm_err]).T
np.savetxt(savepath+archivo_out, dataexport)