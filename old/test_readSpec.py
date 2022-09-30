import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *



# directorio de guradado
savepath  = "S:/PosDoc/Medidas/2022-02-10_gC3N4/"

# directorio de datos
path  = "S:/PosDoc/Medidas/2022-02-10_gC3N4/"

# expnum = [40,42,44,46]
# filename = ['CMK3_1H_5kHz','CMK3_1H_2kHz','CMK3_1H_0kHz','SAC_1H_0kHz']
# expnum = [41,43,45,47]
# filename = ['CMK3_7Li_5kHz','CMK3_7Li_2kHz','CMK3_7Li_0kHz','SAC_7Li_0kHz']
expnum = [108,110]
filename = ['Li2S6_gC3N4', 'Li2S6' ]


i=0
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS
    spec = datos.espectro.spec  
    ppmAxis = datos.espectro.ppmAxis
    # recorto
    # spec = spec[ppmAxis>-20]
    # ppmAxis = ppmAxis[ppmAxis>-20]
    # spec = spec[ppmAxis<20]
    # ppmAxis = ppmAxis[ppmAxis<20]
    # # quito linea de base:
    # base = (spec[0]+spec[-1])/2
    # spec = spec-base        
    
    # vf = VoigtFit(ppmAxis,spec, Npicos=6)
    # ajuste, componentes = vf.componentes(ppmAxis)
    # integrales.append(integrate.trapz(spec))
    # integrales_fit.append(integrate.trapz(ajuste))
    

    dataexport = np.array([ppmAxis, np.real(spec), np.imag(spec)]).T
    np.savetxt('{}{}.dat'.format(savepath,filename[i]), dataexport)
    i+=1
    
    print('graficando...')
    # plt.figure(expn)
    plt.plot(ppmAxis, spec)
    # plt.plot(ppmAxis, ajuste, 'k')
    # for comp in componentes:
    #   plt.plot(ppmAxis, comp, '--')
    # plt.xlim([10,-10])
    # plt.ylim([0,400])

#%%   
