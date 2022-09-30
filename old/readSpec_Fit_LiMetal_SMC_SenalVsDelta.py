"""
Estudio la senal del litio metalico en funcion de delta de la secuencia SMC 16

Estimo barra de error en la determinacion de la amplitud
"""


import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *


# SMC test
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/"
expnum  = [30,31,32,33]
deltas = [2,5,10,20]
filenames = ["SMC16_k1_Delta02ms","SMC16_k1_Delta05ms","SMC16_k1_Delta10ms","SMC16_k1_Delta20ms",]
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/vsDelta/"
# savepath = "S:/temp/"

N = 100
ppm_ini = np.random.uniform(low=100, high=200, size=(N))
ppm_fin = np.random.uniform(low=300, high=400, size=(N))


integrales = np.zeros([len(expnum), 4]).astype('complex')
jj=0
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)
        
    
    real = datos.espectro.real    
    imag = datos.espectro.imag
    spec = real + 1j* imag
    ppmAxis = datos.espectro.ppmAxis
    p = ppmAxis
    # integro en distintos rangos:

    integrales_expn = []
    for nn in range(N):            
      ppm_i = ppm_ini[nn]
      ppm_f = ppm_fin[nn]
      ppm_med = (ppm_f+ppm_i)/2
      ppm_dif = (ppm_f-ppm_i)/2    
      ## recorto a la region deseada
      condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
      spec = spec[condicion]       
      imag = imag[condicion]       
      ppmAxis = ppmAxis[condicion]    
      
      integral =  integrate.simps((spec), x=-ppmAxis) # ABS
      integrales_expn.append(integral)
      
    integrales_expn = np.array(integrales_expn)
    plt.figure(jj)
    plt.hist(np.imag(integrales_expn))
    
    media = np.mean(integrales_expn)
    std   = np.std(np.real(integrales_expn)) + 1j* np.std(np.imag(integrales_expn))
    maxx  = np.max(np.real(integrales_expn)) + 1j* np.max(np.imag(integrales_expn))
    minn  = np.min(np.real(integrales_expn)) + 1j* np.min(np.imag(integrales_expn))
    
    integrales[jj,:] = [media, std, maxx, minn]
        
    
    # # print('graficando...')    
    plt.figure(6578329232)
    plt.plot(ppmAxis, spec)
    plt.plot(ppmAxis, imag,'--')
    

    
    jj+=1


#%%   

media, std, maxx, minn = integrales.T


re  = np.real(media)
im  = np.imag(media)
amp = np.abs(media)
amp_max = np.max(amp)
# normalizo
ren  = re/amp_max
imn  = im/amp_max
ampn = amp/amp_max

std_r = np.real(std)
std_i = np.imag(std)
std_a = np.abs(std)

err_r = ren  * np.sqrt( (std_r/re )**2 + (std_a[amp==amp_max]/amp_max)**2)
err_i = imn  * np.sqrt( (std_i/im )**2 + (std_a[amp==amp_max]/amp_max)**2)
err_a = ampn * np.sqrt( (std_a/amp)**2 + (std_a[amp==amp_max]/amp_max)**2)

err_r = np.max(err_r)
err_i = np.max(err_i)

err_a = np.max(err_a)

plt.figure(123333)
# plt.errorbar(deltas, media_n, yerr=error, marker= 'ko-')
plt.errorbar(deltas, ampn, yerr=err_a, marker= 'o', color='k')
plt.errorbar(deltas, ren, yerr=err_r, marker= 'o', color='b')
plt.errorbar(deltas, imn , yerr=err_i, marker= 'o', color='r' )

# plt.plot(deltas, media_n,  'ko-')
# plt.plot(deltas, re/np.abs(media[0]),  'ro-')
# plt.plot(deltas, im/np.abs(media[0]),  'bo-')

# plt.plot(deltas, amplitudes/amplitudes[0], 'o-')

