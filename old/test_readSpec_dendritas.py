import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
# from VoigtFit import *

# directorio de datos
path = "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/"
path = "S:/Doctorado/LiMetal/116MHz/2021-07-16_litio_in-situ/"
path = "S:/Doctorado/LiMetal/116MHz/2021-10-06_litio_in-operando/"





expnum  = [1,2,3,4]
muestra = [str(exp) for exp in expnum]

expnum  = [1,5,6,7,8,9,10 ]
tiempos = [0, 10,20,30,45,60]
muestra = ["{}min".format(i) for i in tiempos]
muestra.append("60min_NS2048")

expnum  = [5]
tiempos = [10]
muestra = ["{}min".format(i) for i in tiempos]

# cellda T3
path = "S:/Doctorado/LiMetal/116MHz/2021-08-13_litio_in-situ/"
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2021-10_T3_Fityk_BatchFitting/datos2/"
expnum  = np.arange(1,16)
tiempos = [0,1,4,7,13,21,31,45,60,75,90,120,150,180,210]
muestra = ["{}min".format(i) for i in tiempos]



# selector de rango
ppm_i = 100
ppm_f = 400
ppm_med = (ppm_f+ppm_i)/2
ppm_dif = (ppm_f-ppm_i)/2

# inicio un contador
cnt = 0
for expn in expnum:    
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    
    spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    p = ppmAxis
    # recorto a la region deseada
    condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
    spec = spec[condicion]        
    ppmAxis = ppmAxis[condicion]    
    
    archivo_out = str(15-cnt) + "_" + muestra[cnt]
    dataexport = np.array([ppmAxis, spec]).T
    np.savetxt(savepath+archivo_out+'.dat', dataexport)
    
    print('graficando...')
    plt.figure(0)
    plt.plot(ppmAxis, spec)
    
    cnt += 1


ax = plt.gca()
ax.invert_xaxis()
plt.show()