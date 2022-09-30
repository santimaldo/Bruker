import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

import lmfit as lm
import scipy.integrate as integrate

#-------------------------- EXTRAER FACTOR VD----------------------------------
def extraer_vd(path):
    pulseprog = path + 'pulseprogram'
    data = []
    with open(pulseprog, 'rt') as f:    
        for line in f:
          if line.lstrip().startswith('vd*'):
            factor_vd = line.rstrip().split('*')
            print('factor_vd = '+ factor_vd[1])
            factor_vd = float(factor_vd[1])
          else:
            continue
    return factor_vd
#------------------------------------------------------------------------------


# directorio de guradado
path1 = "S:/Doctorado/Carbones/analisis/2019-08_Liberacion/CM7/Remedicion/"
archivo_out = 't_24h'

expn = 18
rser = -1

path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"+str(expn)+"/"
path  = "S:/Doctorado/Carbones/300MHz/2019-09-26_Carbones_Liberacion_CM7/"+str(expn)+"/"


datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real[rser]/datos.acqus.NS

# --------- fit
model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model = model_1 + model_2

params_1 = model_1.make_params(amplitude=7e6, center=-6, sigma=10)
params_2 = model_2.make_params(amplitude=7e5, center=-6, sigma=10)
params = params_1.update(params_2)

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

#path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"+str(expn)+"/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real/datos.acqus.NS
maximo = np.max(spec)

model_1 = lm.models.VoigtModel(prefix='m1_')
model_2 = lm.models.VoigtModel(prefix='m2_')
model = model_1 + model_2 

#plt.figure(0)
integrales = np.zeros([len(componentes)+1, len(spec)])
for j in range(len(spec)-1, -1,-1):
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
    for i in range(len(componentes)):                
        comp = componentes[f'm{i+1}_']
        plt.plot(ppmAxis, comp, '--')
        integrales[i, j] = integrate.trapz(comp, ppmAxis)
    integrales[i+1, j] = integrate.trapz(spec[j], ppmAxis)
#%%
tau = np.loadtxt(path+'vdlist')*1000 # paso a milisegundos
npts = integrales[0].size
tau=tau[0:npts]
tau = tau * extraer_vd(path)

S1 = -integrales[0]
S2 = -integrales[1]
S = S1 + S2
integral_total = -integrales [2]
intensidad = np.max(spec, axis=1)
    
plt.figure(1)
plt.plot(tau, S1, 'bo-')
plt.plot(tau, S2, 'ro-')
plt.plot(tau, S, 'ko-')
plt.plot(tau, intensidad, 'k--')
plt.plot(tau, integral_total, 'k-')

T1data = np.array([tau, S1, S2, S, intensidad]).T

np.savetxt(path1+archivo_out+'.dat', T1data)


plt.show()


