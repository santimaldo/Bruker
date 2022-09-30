import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *



# directorio de guradado
savepath = "S:/Doctorado/Carbones/analisis/2022-Glyme_AnalisisDatos2018/Datos/"



# directorio de datos
# path = "S:/Doctorado/Carbones/300MHz/2018-09-25_Carbones_Glyma/" # Bulks
# path = "S:/Doctorado/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/" # Confinados

# CMK3ACT
savepath = "S:/Doctorado/Carbones/analisis/2021-12_Carbones_MAS/files/spec1D"
path = "S:/Doctorado/Carbones/300MHz/2022-05-06_Carbones_CMK3ACT/" # CMK3ACT
archivos = ['TEGDME_Bulk_1H','TEGDME_Bulk_7Li','CMK3ACT_1H','CMK3ACT_7Li']
expnum = [4,5,13,14]
ppmLim = 20

# CMK3ACT-B
savepath = "S:/Doctorado/Carbones/analisis/2021-12_Carbones_MAS/files/"
path = "S:/Doctorado/Carbones/300MHz/2022-05-12_Carbones_CMK3ACT/" # CMK3ACT
archivos = ['CMK3ACT-B_1H', 'CMK3ACT-B_7Li']
expnum = [14,20]
ppmLim = 30

# CMK3ACT-B
savepath = "S:/temp/"
path = "S:/Doctorado/Carbones/300MHz/2022-05-19_Carbones_CMK3ACT/" # CMK3ACT
# expnum = [1] + [i for i in range(103,107)] + [i for i in range(108,125)] + [1000]
expnum = [1] + [100] + [i for i in range(102,107)] + [i for i in range(108,125)] + [1000]
archivos = []
for i in range(len(expnum)):
  archivos.append(f"1H_nexp{expnum[i]}")

ppmLim = 30


# archivos = ['7Li_LiTF-in-diglyme_HF','7Li_LiTF-in-diglyme_NaOH']
# expnum = [4,10]

# archivos = ['1H_LiTF-in-diglyme_HF','1H_LiTF-in-diglyme_NaOH']
# expnum = [3,9]
# ppmLim = None

rangos = [[-4,-3], [2,2.6], [3.4,3.8]]
Integrales = []

Tiempo = []

Amp_list = []
nn = -1
for expn in expnum:    
    nn+=1
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    Tiempo.append(datos.acqus.hora)
    
    spec = datos.espectro.real / NS    
    ppmAxis = datos.espectro.ppmAxis
    # recorto el eje entre -ppmLim y ppmLim
    if ppmLim is not None:
      spec = spec[np.abs(ppmAxis)<ppmLim]
      ppmAxis = ppmAxis[np.abs(ppmAxis)<ppmLim]
        
    
    ## integro:
    for n_rango in range(len(rangos)):
      rango = rangos[n_rango]
      # defino arrays:
      if nn==0:
        Integrales.append(np.zeros_like(np.array(expnum)))
      
      ppmi=rango[0]; ppmf=rango[1]
      
      spec_tmp = spec[np.abs(ppmAxis-0.5*(ppmi+ppmf))<0.5*np.abs(ppmf-ppmi)]
      ppmAxis_tmp = ppmAxis[np.abs(ppmAxis-0.5*(ppmi+ppmf))<0.5*np.abs(ppmf-ppmi)]
      Integrales[n_rango][nn] = np.abs(integrate.simps(spec_tmp, x=ppmAxis_tmp))
    
    amplitud = np.abs(integrate.simps(spec, x=ppmAxis))
    Amp_list.append(amplitud)
    
    archivo_out = archivos[nn]
    dataexport = np.array([ppmAxis, spec, spec/amplitud]).T
    np.savetxt(savepath+archivo_out+'.dat', dataexport)
    
    print('graficando...')
    plt.figure(111110)
    plt.plot(ppmAxis, spec/amplitud)



plt.plot(ppmAxis, spec/amplitud*0, 'k--')
#%% 
# t = [0,1]+[0.5*i+1 for i in range(1)] + [10+i*2 for i in range(2)] + [12+i*4 for i in range(4)] + [34,44]
t = np.array(Tiempo)
t = t - t[0]
t = t/60 # paso a minutos
t[-1] = 300 # timpo arbitrario pero grande para el punto "final"
n=0
for integral in Integrales:
  plt.figure(12555)
  plt.plot(t, integral,'o-', label=f'{rangos[n]} ppm')
  n+=1      
plt.legend()
plt.xlabel("Tiempo [min]")

Amp = np.array(Amp_list)
plt.figure(231)
plt.plot(t, Amp/np.max(Amp), 'o-')
plt.xlabel("Tiempo [min]")
# # espectro sobre el cual normalizar        
# directorio = path+str(1)+"/"
# datos = DatosProcesados(directorio)    
# NS = datos.acqus.NS  

# spec = datos.espectro.real / NS    
# ppmAxis = datos.espectro.ppmAxis
# #    ppmAxis = ppmAxis[ppmAxis>100]
# #    re = re[0:ppmAxis.size]        
# vf = VoigtFit(ppmAxis,spec, Npicos=2)
# ajuste, componentes = vf.componentes(ppmAxis)

# FactorNorm = integrate.trapz(spec)
# FactorNorm_fit = integrate.trapz(ajuste)

        
# plt.figure(789)
# integrales = integrales/FactorNorm
# integrales_fit = integrales_fit/FactorNorm_fit
# plt.plot(tiempos, integrales,'o-', label='Integral de los datos')
# plt.plot(tiempos, integrales_fit,'o-', label='Integral del ajuste')
# plt.show()