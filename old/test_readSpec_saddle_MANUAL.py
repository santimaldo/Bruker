import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
# from VoigtFit import *

# directorio de datos
path = "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/"
# directorio de guradado
savepath = "S:/temp/"



expnum = [3*i+1 for i in range(9)]
muestra = [str(exp) for exp in expnum]

rho = np.array([0,6,6,3,3,2,4,6,3])
phi = np.array([0,180,90,180,90,90,90,90,180])
phi_r = phi*np.pi/180
#------------------------------------------------------------
expnum = [61,121,124]
rho = np.array([0,0,0])
phi = np.array([0,0,0])
phi_r = phi*np.pi/180

Amp = []
Ph = []
# inicio un contador
cnt = 0
for expn in expnum:   
  directorio = path+str(expn)+"/"
  datos = DatosProcesados(directorio)    
  NS = datos.acqus.NS  
  ph = datos.procs.phase
  
  # spec = datos.espectro.real + 1j * datos.espectro.imag
  spec = datos.espectro.real
  ppmAxis = datos.espectro.ppmAxis
  
  ini = 3
  fin = -5.5
  spec = spec[ppmAxis<ini]
  ppmAxis = ppmAxis[ppmAxis<ini]
  spec = spec[ppmAxis>fin]
  ppmAxis = ppmAxis[ppmAxis>fin]

  # plt.figure(expn)
  # plt.plot(ppmAxis, spec)
  # plt.xlim([20,-11])
  
  # integral = np.trapz(np.abs(spec))
  integral = np.trapz(spec)
  Amp.append(integral)
  Ph.append(ph)
  
  # archivo_out = muestra[cnt]
  # dataexport = np.array([ppmAxis, spec]).T
  # np.savetxt(savepath+archivo_out+'.dat', dataexport)
  
  print('graficando...')
  plt.figure(0)
  plt.plot(ppmAxis, spec)
  plt.xlim([20,-11])
  
  cnt += 1

#%%
datos = np.array([rho, phi, Amp,Ph]).T
print(datos)
np.savetxt(savepath+'amplitudes_manual.dat', datos, fmt='%i')


