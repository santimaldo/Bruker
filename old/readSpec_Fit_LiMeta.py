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

# SMC test
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/"
expnum  = [105, 106, 107, 111,112,113,114]
filenames = ["SP_NS64", "SMC16_k1_NS128", "SMC16_k2_NS256", "SMC16_k1","SMC32_k1","SMC64_k1","SP" ]
# expnum  = [30,31,32,33]
# filenames = ["SMC16_k1_Delta02ms","SMC16_k1_Delta05ms","SMC16_k1_Delta10ms","SMC16_k1_Delta20ms"]
expnum  = [51,52,53]
filenames = ["SMC16_k1", "SP", "SMC16_k2"]
expnum  = [114]
filenames = ["SP" ]
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/"


# Dendritas vs Angulo
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2019-04-05_Dendritas/"
expnum   = np.array([1]+[i for i in np.arange(2,48,2)])
filenames = [f"LiMetal_Dend_{i}.dat" for i in expnum]
# directorio de guradado
savepath  = "S:/Doctorado/LiMetal/Analisis/2019-04_TestPortamuestra/Espectro_vs_t-LiDendritas/"
# # Li Metal vs Angulo
# # directorio de datos
# path  = "S:/Doctorado/LiMetal/116MHz/2019-04-17_Dendritas/"
# expnum   = [i for i in np.arange(2,44,2)]
# filenames = [f"LiMetal_{i}.dat" for i in expnum]
# # directorio de guradado
# savepath  = "S:/Doctorado/LiMetal/Analisis/2019-04_TestPortamuestra/Espectro_vs_t/"


# Dendritas vs Angulo
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2019-04-16_Dendritas/"
expnum   = [8]
filenames = [f"8"]
# directorio de guradado
savepath  = path


# selector de rango
ppm_i = -100
ppm_f = 400
ppm_med = (ppm_f+ppm_i)/2
ppm_dif = (ppm_f-ppm_i)/2

plt.figure(123568)
plt.xlim((ppm_f,ppm_i))

Tiempo = []

integrales = []
integrales_err = []
jj=0
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)
    
    spec = datos.espectro.real
    imag = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis
    p = ppmAxis
    Tiempo.append(datos.acqus.hora)
    # recorto a la region deseada
    condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
    spec = spec[condicion]       
    imag = imag[condicion]       
    ppmAxis = ppmAxis[condicion]    

    integrales.append(np.abs(integrate.simps(spec+1j*imag, x=ppmAxis))) # ABS
    # integrales.append(np.abs(integrate.simps(imag, x=ppmAxis))) # 

    archivo_out = filenames[jj]+'.dat'
    dataexport = np.array([ppmAxis, spec, imag, spec/np.max(spec), imag/np.max(spec)]).T
    np.savetxt(savepath+archivo_out, dataexport)
    # print('graficando...')    
    plt.plot(ppmAxis, spec)
    # plt.plot(ppmAxis, imag,'--')
    
    jj+=1

integrales = np.array(integrales)    
integrales_err = np.array(integrales_err)
#%%   
t = np.array(Tiempo)
t = t - t[0]
t = t/60/60 # paso a horas

plt.figure(1)
plt.plot(integrales)
plt.xlabel("Tiempo [h]")