import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

import lmfit as lm
import scipy.integrate as integrate

expn = 2000
path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_SiO2/"+str(expn)+"/"

expn = 5
path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"+str(expn)+"/"



datos = DatosProcesados(path)

spec = datos.espectro.real
ppmAxis = datos.espectro.ppmAxis
#ppmAxis = ppmAxis[ppmAxis>100]
#spec = spec[0:ppmAxis.size]


# --------- fit
#model_1 = lm.models.VoigtModel(prefix='m1_')
#model_2 = lm.models.VoigtModel(prefix='m2_')
#model = model_1 + model_2

Npicos = 2
model = lm.models.VoigtModel(prefix='m1_')
params = model.make_params(amplitude=1e8, center=245, sigma=1, gamma=1)

models_i = [model]
params_i = [params]
for i in range(1,Npicos):
    prefijo_i = 'm{}_'.format(i+1)
    model = model + lm.models.VoigtModel(prefix=prefijo_i)
    params = model.make_params(amplitude=1e8, center=261, sigma=1, gamma=1)
    models_i.append(model)
    params_i.append(params)
    
#
#params_1 = models_i[0].make_params(amplitude=1e8, center=245, sigma=10)
#params_2 = models_i[1].make_params(amplitude=1e7, center=260, sigma=10)
params = params_i[0].update(params_i[1])

ajuste = model.fit(spec, params, x=ppmAxis)
componentes = ajuste.eval_components(x=ppmAxis)

plt.figure(456)
plt.plot(ppmAxis, spec, 'b')    
plt.plot(ppmAxis, ajuste.best_fit, 'k-', label='best fit')
for i in range(1,len(componentes)+1):
    plt.plot(ppmAxis, componentes[f'm{i}_'], '--')

plt.legend(loc='best')
plt.show()
