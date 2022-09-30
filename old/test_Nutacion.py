import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

import lmfit as lm
import scipy.integrate as integrate



path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_Nutacion/2/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = (datos.espectro.real + 1j*datos.espectro.imag)/datos.acqus.NS

ppmRange = [244,246]
ppmAxis = datos.espectro.ppmAxis
ppmAxis = ppmAxis[ppmAxis>ppmRange[0]]
spec = spec[:,0:ppmAxis.size]
ppmAxis = ppmAxis[ppmAxis<ppmRange[1]]
spec = spec[:,-1:-ppmAxis.size-1:-1]
 

real = np.real(spec)
imag = np.imag(spec)


#plt.figure(0)
integrales_re = []
integrales_im = []
plt.figure(0)
for j in range(len(real)):    
    integral_re = integrate.trapz(real[j])
    integral_im = integrate.trapz(imag[j])
    integrales_re.append(integral_re)    
    integrales_im.append(integral_im)
    plt.plot(ppmAxis, real[j])

#cuantos pasos en tp hubo:
tp1 = np.linspace(1,len(real),len(real))            
#%%
            
path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_Nutacion/4/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = (datos.espectro.real + 1j*datos.espectro.imag)/datos.acqus.NS

ppmAxis = datos.espectro.ppmAxis
ppmAxis = ppmAxis[ppmAxis>ppmRange[0]]
spec = spec[:,0:ppmAxis.size]
ppmAxis = ppmAxis[ppmAxis<ppmRange[1]]
spec = spec[:,-1:-ppmAxis.size-1:-1]
 

real = np.real(spec)
imag = np.imag(spec)


#plt.figure(0)
plt.figure(0)
for j in range(len(real)):    
    integral_re = integrate.trapz(real[j])
    integral_im = integrate.trapz(imag[j])
    integrales_re.append(integral_re)    
    integrales_im.append(integral_im)
    plt.plot(ppmAxis, real[j])            

#%%
tp2 = np.linspace(128,128+31,32)
tp = np.concatenate((tp1, tp2))

    
plt.figure(1)
plt.plot(tp, integrales_re, 'ko-')
plt.plot(tp, integrales_im, 'ro-')

#T1data = np.array([tp, S1, S2, S, intensidad]).T

#np.savetxt(path1+'13h.dat', T1data)


plt.show()


