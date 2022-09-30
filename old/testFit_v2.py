import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

import lmfit as lm
import scipy.integrate as integrate

# directorio de guradado
path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_Nutacion/2/"


datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real[15]/datos.acqus.NS

ppmAxis = datos.espectro.ppmAxis
#ppmAxis = ppmAxis[ppmAxis>10]
#spec = spec[0:ppmAxis.size]

# --------- fit
model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model_3 = lm.models.VoigtModel(prefix='m3_')
model = model_1 + model_2 + model_3

params_1 = model_1.make_params(amplitude=7e6, center=-6, sigma=10)
params_2 = model_2.make_params(amplitude=7e5, center=-6, sigma=10)
params_3 = model_3.make_params(amplitude=7e5, center=-6, sigma=10)
params = params_1.update(params_2)
params = params_1.update(params_3)

ajuste = model.fit(spec, params, x=ppmAxis)
componentes = ajuste.eval_components(x=ppmAxis)

plt.figure(456)
plt.plot(ppmAxis, spec, 'b')    
plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')
for i in range(1,len(componentes)+1):
    plt.plot(ppmAxis, componentes[f'm{i}_'], '--')

plt.legend(loc='best')

plt.figure(0)
plt.subplot(6,6,36)
plt.plot(ppmAxis, spec, 'b')    
plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')
for i in range(1,len(componentes)+1):
    plt.plot(ppmAxis, componentes[f'm{i}_'], '--')

plt.legend(loc='best')




#%%
#---------------------- 2D

path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_Nutacion/2/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real/datos.acqus.NS
maximo = np.max(spec)
ppmAxis = datos.espectro.ppmAxis
ppmAxis = ppmAxis[ppmAxis>100]
spec = spec[0:ppmAxis.size]


model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model_3 =  lm.models.VoigtModel(prefix='m3_')
model = model_1 + model_2 + model_3

#plt.figure(0)
integrales = np.zeros([len(componentes), len(spec)])
for j in range(len(spec)):
    params = ajuste.params
    for i in range(1,len(componentes)+1):
        params[f'm{i}_center'].vary = False
        params[f'm{i}_amplitude'].min = 0
        
    
    ajuste = model.fit(spec[j], params, x=ppmAxis)
    componentes = ajuste.eval_components(x=ppmAxis)        
    plt.subplot(6, 6,j+1)    
    plt.plot(ppmAxis, spec[j], 'b')    
    plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')    
    plt.ylim([-0.1*maximo, maximo])
    plt.legend(loc='best')
    if (i==15 or i==30 or i==44):    
        for i in range(len(componentes)):                
            comp = componentes[f'm{i+1}_']
            plt.plot(ppmAxis, comp, '--')
            integrales[i, j] = integrate.trapz(comp, ppmAxis)             
#%%
tp = np.linspace(1,128,128)            
S1 = -integrales[0]
S2 = -integrales[1]
S = S1 + S2 
intensidad = np.max(spec, axis=1)
    
plt.figure(1)
plt.plot(tp, S1, 'bo-')
plt.plot(tp, S2, 'ro-')
plt.plot(tp, S, 'ko-')
plt.plot(tp, intensidad, 'k--')

#T1data = np.array([tp, S1, S2, S, intensidad]).T

#np.savetxt(path1+'13h.dat', T1data)


plt.show()


