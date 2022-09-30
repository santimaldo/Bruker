import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

import lmfit as lm
import scipy.integrate as integrate

path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/13/"
path = "S:/Doctorado/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/5/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

# --------- fit
spec=spec[-1]

model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model_3 = lm.models.VoigtModel(prefix='m3_')
model_4 = lm.models.VoigtModel(prefix='m4_')
model = model_1 + model_2 + model_3 + model_4

params_1 = model_1.make_params(amplitude=7e6, center= 0, sigma=4)
params_2 = model_2.make_params(amplitude=7e5, center= 0, sigma=4)
params_3 = model_3.make_params(amplitude=7e5, center=-6, sigma=5)
params_4 = model_4.make_params(amplitude=7e5, center= 0 , sigma=50)
params = params_1.update(params_2)
params = params_1.update(params_3)
params = params_1.update(params_4)

ajuste = model.fit(spec, params, x=ppmAxis)
componentes = ajuste.eval_components(x=ppmAxis)

plt.figure(456)
plt.plot(ppmAxis, spec, 'b')    
plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')
for i in range(1,len(componentes)+1):
    plt.plot(ppmAxis, componentes[f'm{i}_'], '--')

plt.legend(loc='best')
plt.show()

#%%
#---------------------- 2D

path = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/14/"
path = "S:/Doctorado/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/5/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model_3 = lm.models.VoigtModel(prefix='m3_')
model = model_1 + model_2 + model_3

plt.figure(0)
integrales = np.zeros([3, 32])
jj=0
for j in range(31,27,-1):
    params = ajuste.params
    for i in range(1,len(componentes)+1):
        params[f'm{i}_center'].vary = False
        params[f'm{i}_amplitude'].min = 0
        
    
    ajuste = model.fit(spec[j], params, x=ppmAxis)
    componentes = ajuste.eval_components(x=ppmAxis)
    
    jj += 1
    plt.subplot(1, 4,jj)    
    plt.plot(ppmAxis, spec[j], 'b')    
    plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')
    plt.xlim([-20, 20])
    plt.ylim([-1e4, 2e5])
    
    
    for i in range(len(componentes)):                
        comp = componentes[f'm{i+1}_']
        plt.plot(ppmAxis, comp, '--')
        integrales[i, j] = integrate.trapz(comp, ppmAxis)        
    
    plt.legend(loc='best')
    
#%%
tau = np.loadtxt(path+'vdlist')
    
plt.figure(1)
plt.plot(tau, -integrales[0], 'bo-')
plt.plot(tau, -integrales[1], 'ro-')
plt.plot(tau, -integrales[2], 'go-')
plt.plot(tau, -integrales[1]-integrales[0]-integrales[2], 'ko--')


    
plt.show()


