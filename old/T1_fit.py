import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

from VoigtFit import *
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
path1 = "S:/tmp/"
archivo_out = 'tmp'

expn = 2
rser = -1

path = f"S:/Carbones/300MHz/2022-09-01_CarbonesCNEA_M4_LiTFSI/{expn}/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real[rser]/datos.acqus.NS

# --------- fit

vf = VoigtFit(ppmAxis, spec, Npicos=2, center=[0,0])
ajuste, componentes = vf.componentes(ppmAxis)
parametros_ini = vf.params

plt.figure(456)
plt.plot(ppmAxis, spec, 'b')
plt.plot(ppmAxis, ajuste, 'k-', label='best fit')
for i in range(len(componentes)):
    plt.plot(ppmAxis, componentes[i], '--')

plt.legend(loc='best')

plt.figure(0)
plt.subplot(6,6,36)
plt.plot(ppmAxis, spec, 'b')
plt.plot(ppmAxis, ajuste, 'k-', label='best fit')
for i in range(len(componentes)):
    plt.plot(ppmAxis, componentes[i], '--')

plt.legend(loc='best')




#%%
#---------------------- 2D

#path  = "S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/"+str(expn)+"/"

datos = DatosProcesados(path)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real/datos.acqus.NS
maximo = np.max(spec)
vf = None
#plt.figure(0)
integrales = np.zeros([len(componentes)+1, len(spec)])
for j in range(len(spec)-1, -1,-1):
    vf = VoigtFit(ppmAxis, spec, Npicos=2, params=parametros_ini, fijar= 'center')
    #vf = VoigtFit(ppmAxis, spec, Npicos=2, params=params)
    ajuste, componentes = vf.componentes(ppmAxis)

    plt.subplot(6, 6,j+1)
    plt.plot(ppmAxis, spec[j], 'b')
    plt.plot(ppmAxis, ajuste, 'k-', label='best fit')
    plt.ylim([-0.1*maximo, maximo])
    plt.legend(loc='best')
    for i in range(len(componentes)):
        comp = componentes[i]
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


