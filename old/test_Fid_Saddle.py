import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Espectro import *

path  = path = "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/"
savepath = "S:/temp/"


r =      [-6 ,-5  ,-4 ,-3 ,-2 ,-1 ,0  ,1 ,2 ,3 ,4 ,5  ,6 ]
expnum = [ 82, 106, 52, 58, 49, 94,121,88,40,34,37,100,31] #x
expnum = [ 43, 112, 73, 46, 67, 97,121,91,70,85,76,103,79] #y


integral=[]
fase=[]

for expn in expnum:

  datos = DatosCrudos(path+str(expn))
  
  fid = datos.fid
  ###adquiero menos tiempo: menos ruido
  # tmax = 20e-3
  # fid.RecortarTiempo(tmax)
  # # # exponential multiplication:
  fid.em(100)
  
  espectro = Espectro()
  
  ppmAxis, spec, angle = espectro.CrearEspectro(fid, return_angle=True)
  
  spec = spec/datos.acqus.NS/datos.acqus.RG
  mod = np.abs(spec)
  mod = mod - np.mean(mod[300:])
  
  plt.figure(1)
  # plt.plot(ppmAxis, mod)
  plt.plot(mod)
  plt.plot(mod*0,'k--')
  
  integral.append(np.trapz(mod))
  fase.append(angle)

#%%
integral = [i/integral[6] for i in integral]
plt.figure(0)
plt.plot(r,integral,'o-')
plt.figure(10)
ph = [f-fase[6] for f in fase]
plt.plot(r,ph,'o-')