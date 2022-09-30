import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
# from VoigtFit import *

# directorio de datos
# path = "S:/Doctorado/LiMetal/116MHz/2021-11-04_in-situ_T10/"
# path = "S:/Doctorado/LiMetal/116MHz/2021-11-10_in-situ_T11-12/" # T11 a T13
# directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_1carga_T9_a_T13/"
# experimentos:
expnum  = np.arange(1,16)
# nombres de archivo
# T9 -------------------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2021-11-03_in-situ_T9/" 
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_1carga_T9_a_T13/T9/"
expnum  = [1,2,3 ,4 ,105,7 ,8 ]
tiempos = [0,5,10,20,35 ,50,60]
archivo = ["T9_{}min".format(i) for i in tiempos]

# T10 ------------------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2021-11-04_in-situ_T10/" 
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_1carga_T9_a_T13/T10/"
expnum  = [1, 2]
tiempos = [0,60]
archivo = ["T10_{}min".format(i) for i in tiempos]

# T11 a 13------------------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2021-11-10_in-situ_T11-12/" 
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_1carga_T9_a_T13/"
expnum  = [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 11]
archivo = ["T11_0min_A","T11_0min_B","T11_60min_A","T11_60min_B",\
            "T12_0min_A","T12_0min_B","T12_60min_A","T12_60min_B",\
            "T13_0min","T13_60min_A","T13_60min_B"]
# caract saddle----------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2021-11-25_Homogeneidad_saddle/"  
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_CaracterizacionSaddle/datos/"
expnum  = np.arange(46,66)
archivo = ["{}us_C","{}us_B","{}us_D","{}us_A","{}us_E"]*4


# SMC Li Pristino-------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/"  
savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/LiPristino/"
expnum  = np.arange(51,54)
archivo = ["SMC_k1", "SP", "SMC_k2"]

# SMC Li Ciclado-------------------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/"  
savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/LiCiclado/"
expnum  = [114,106,107, 111, 112, 113]
archivo = ["SP", "SMC_k1",  "SMC_k2", "SMC_k1_N16", "SMC_k1_N32", "SMC_k1_N64"]

# SMC Li Pristino vs Delta-----------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/"  
savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/vsDelta/"
expnum  = [30,31,32,33]
archivo = ["SMC16_k1_Delta02ms", "SMC16_k1_Delta05ms", "SMC16_k1_Delta10ms", "SMC16_k1_Delta20ms"]

# SiO2-----------------------------------------------------
path = "S:/Doctorado/LiMetal/116MHz/2020-03-12_Dendritas_SiO2/"  
savepath = "S:/Doctorado/LiMetal/Analisis/2019-09_SiO2/"
expnum  = [1,2,3,4]
archivo = ["M_0.5mAcm2_ctrl_E+", "M_0.5mAcm2_ctrl_E-", "M_0.5mAcm2_SiO2_E-", "M_0.5mAcm2_SiO2_E+"]


Integrales = True

# selector de rango
ppm_i = 100
ppm_f = 400
ppm_med = (ppm_f+ppm_i)/2
ppm_dif = (ppm_f-ppm_i)/2

Areas = []
tp = []
# inicio un contador
cnt = 0
for expn in expnum:    
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)    
    
    spec = datos.espectro.real
    imag = datos.espectro.imag
    
    ppmAxis = datos.espectro.ppmAxis
    P1 = datos.acqus.P1

    # recorto a la region deseada

    condicion = np.abs(ppmAxis-ppm_med)<=np.abs(ppm_dif)
    spec = spec[condicion]        
    imag = imag[condicion]        
    ppmAxis = ppmAxis[condicion]    
    
    
    # dataexport = np.array([ppmAxis, spec, imag, spec/np.max(spec), imag/np.max(spec)]).T
    # dataexport = np.array([ppmAxis-7, spec, spec/np.max(spec)]).T
    dataexport = np.array([ppmAxis, spec, imag]).T
    np.savetxt(f"{savepath}{archivo[cnt]}.dat", dataexport)
    
    if Integrales:
      area = np.abs(integrate.simps(spec, x=ppmAxis))
      Areas.append(area)
      tp.append(P1)
    print('graficando...')
    plt.figure(0)
    plt.plot(ppmAxis, spec)
    plt.plot(ppmAxis, imag)
      
    
    cnt += 1


plt.xlim([ppm_f,ppm_i])
plt.xlabel("NMR shift [ppm]")
plt.show()

#%%
if Integrales:
  Areas = np.array(Areas)
  tp = np.array(tp)
  plt.figure(111)
  plt.plot(Areas, 'o-')

  
  
  
  
  