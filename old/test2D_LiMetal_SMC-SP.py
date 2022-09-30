import nmrglue as ng
import matplotlib.pyplot as plt
from scipy.integrate import simps
import numpy as np
from Datos import *
from Autophase import autophase
plt.rcParams.update({'font.size': 20})
"""
Al SMC vs k le resto el espectru bulk

deben tener ppmAxis identicos
"""


#
pathBulk = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/21/" # bulk 1D

path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1022/"  # pristino
k_list = (np.arange(32) + 4.5)*2 /14
filename = "LiPristino"

# path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1020/" # pristino
# k_list = (np.arange(16)*0.25 + 6)*2 /14
# filename = "LiPristino1"

# path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1110/" # Li Ciclado
# k_list = (np.arange(19)*0.2 + 5.3)*2 /14
# filename = "LiCiclado2"

path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1109/" # Li Ciclado
k_list = (np.arange(64)*0.5 + 4.5)*2 /14
filename = "LiCiclado1"


datos = DatosProcesados(path)
datosBulk = DatosProcesados(pathBulk)

# limites de los specs
ppm_i = 220
ppm_f = 280
ppm_med = (ppm_f+ppm_i)/2
ppm_dif = (ppm_f-ppm_i)/2    

# extraigo spec 2D
ppmAxis = datos.espectro.ppmAxis
npts= datos.espectro.size[0]
## recorto a la region deseada
condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
ppmAxis = ppmAxis[condicion]    
resolucion = ppmAxis[0]-ppmAxis[1]


# Defino el bulk como el k=2
# ind = np.where(k_list==2)[0][0]
ind = np.where(np.abs(k_list-2)<0.1)[0][0]


ppmAxisB = datos.espectro.ppmAxis
realB = datos.espectro.real[ind]
imagB = datos.espectro.imag[ind]    

## recorto a la region deseada
condicion = np.abs(ppmAxisB-ppm_med)<np.abs(ppm_dif)
ppmAxisB = ppmAxisB[condicion]
realB = realB[condicion]
imagB = imagB[condicion]
realB, imagB, phase = autophase(ppmAxisB, realB, imagB)

magB = np.abs(realB+1j*imagB)




#%%

realInt = np.zeros(npts)
imagInt = np.zeros(npts)
substract = np.zeros(npts)
cond = 0
for j in range(npts):    
    real = datos.espectro.real[j]
    imag = datos.espectro.imag[j]    
    ## recorto a la region deseada    
    real = real[condicion]       
    imag = imag[condicion]       
    
    real, imag, phase = autophase(ppmAxis, real, imag)
        
    realInt[j] = simps(real, x=ppmAxis)
    imagInt[j] = simps(imag, x=ppmAxis)
    
    
    # spec = np.abs(real+1j*imag)
    # specB = magB
    
    spec = real
    specB = realB
    
    ## Resto el bulk
    ppmBulk = ppmAxisB[specB==np.max(specB)]
    
    ppmBulk = 252
    factorNorm  = np.mean( spec[np.abs(ppmAxis-ppmBulk)<2]) # valor del espectro en la posicion del bulk
    factorNormB = np.mean(specB[np.abs(ppmAxis-ppmBulk)<2]) # valor del espectro en la posicion del bulk
    
    spec  =  spec/factorNorm
    specB = specB/factorNormB
        
    substract[j] = -simps(spec-specB, x=ppmAxis)
    

    if k_list[j]>0.95:
      if cond == 0:
        jtmp = j
        cond = 1
      if cond == 1:
        if j-jtmp<9:
          plt.figure(789555)
          plt.subplot(3,3,(j-jtmp)+1)
          plt.title(f"k={k_list[j]:.2f}")
          plt.plot(ppmAxis, spec)
          plt.plot(ppmAxis, specB)
          # plt.plot(ppmAxis, specB*factorNorm)
          plt.plot(ppmAxis, spec-specB)
          # plt.xticks([])
          # plt.yticks([])
          plt.xlim([280,230])
    
    
    

#%%


S = realInt + 1j*imagInt

header = 'k    abs(S)    phase(S)    real(S)    imag(S)'
datos = np.array([k_list, np.abs(S), np.angle(S), np.real(S), np.imag(S)]).T
np.savetxt(f"{filename}.dat", datos, header=header)


plt.figure(951)
plt.plot(k_list, np.abs(S)/np.max(np.abs(S)), 'o')
plt.xlabel("k")

# plt.figure(952)
# plt.plot(k_list, np.angle(S)*180/np.pi, 'o')

#%%
plt.figure(953)
plt.plot(k_list[k_list>=1], substract[k_list>=1], 'o-')
plt.xlabel("k")
plt.ylabel(r"$SpecSMC - SpecBulk\times[\frac{SpecSMC(\sim246)}{SpecBulk(\sim246)}]$")
plt.hlines(0,1,5,ls='--', color='k')

plt.show()
