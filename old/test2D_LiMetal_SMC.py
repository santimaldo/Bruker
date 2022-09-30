import nmrglue as ng
import matplotlib.pyplot as plt
from scipy.integrate import simps
import numpy as np
from Datos import *
from Autophase import autophase

# EXTRAIGO BULK
path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1022/"


path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1022/"
k_list = (np.arange(32) + 4.5)*2 /14

# path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1020/"
# k_list = (np.arange(16)*0.25 + 6)*2 /14

savepath = "S:/Doctorado/LiMetal/Analisis/2022-03_SMC_test/BatchFit/EspectrosCIclado/"
path  = "S:/Doctorado/LiMetal/116MHz/2022-03-23_SMC_test/1109/"
k_list = (np.arange(64)*0.5 + 4.5) /14


datos = DatosProcesados(path)


ppmAxis = datos.espectro.ppmAxis
npts= datos.espectro.size[0]


ppm_i = 200
ppm_f = 300
ppm_med = (ppm_f+ppm_i)/2
ppm_dif = (ppm_f-ppm_i)/2    
## recorto a la region deseada
condicion = np.abs(ppmAxis-ppm_med)<np.abs(ppm_dif)
ppmAxis = ppmAxis[condicion]    

realInt = np.zeros(npts)
imagInt = np.zeros(npts)
plt.figure(1)
for j in range(npts):    
    real = datos.espectro.real[j]
    imag = datos.espectro.imag[j]    
    ## recorto a la region deseada    
    real = real[condicion]       
    imag = imag[condicion]       
    
    real, imag, phase = autophase(ppmAxis, real, imag)
    plt.plot(ppmAxis, real)
    
    savedata = np.array([ppmAxis, real, imag]).T
    np.savetxt(f"{savepath}k{k_list[j]:.2f}.dat", savedata)
    
    
    realInt[j] = simps(real, x=ppmAxis)
    imagInt[j] = simps(imag, x=ppmAxis)

#%%
S = realInt + 1j*imagInt

plt.figure(951)
plt.plot(k_list, np.imag(S)/np.max(np.abs(S)), 'ro')
plt.plot(k_list, np.real(S)/np.max(np.abs(S)), 'bo')





plt.show()
