# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:34:06 2022

@author: santi

Corrijo automaticamente la fase, utilizando el método de minimizar el area
de la parte imaginaria:   case{'MinIntImagSpec'}
"""
import numpy as np
from scipy.integrate import simps

def autophase(ppmAxis, real, imag, precision=1):
  
  spec = real + 1j*imag  
  spec = spec[ppmAxis.argsort()] # ordeno basandome en ppmAxis
  ppmAxis = ppmAxis[ppmAxis.argsort()] # ordeno basandome en ppmAxis  
  
  angle=np.arange(-180,180,precision)

  IntImagSpec = []
  IntRealSpec = []
  for i in range(angle.size):
      Sp_try= spec*np.exp(-1j*angle[i]*np.pi/180)      
      IntImagSpec.append(simps(np.imag(Sp_try),x=ppmAxis))
      IntRealSpec.append(simps(np.real(Sp_try),x=ppmAxis))
      
  IntImagSpec = np.array(IntImagSpec)
  IntRealSpec = np.array(IntRealSpec)
  # Ordeno basandome en el modulo de integral imaginaria.
  # luego, si la integral real es negativa, ese angulo no me sirve y paso
  # al siguiente
  dataArray = np.array([np.abs(IntImagSpec), IntRealSpec, angle]).T
  dataArray = dataArray[dataArray[:, 0].argsort()] # ordeno basandome en la primera columna
  
  for i in range(IntImagSpec.size):
    integralReal = dataArray[i,1]
    if integralReal<0:
      # print("NO  ", dataArray[i,:])
      continue
    else:
      # print("SI  ", dataArray[i,:])
      anguloOptimo = dataArray[i,2]
      break    
  spec = spec*np.exp(-1j*anguloOptimo*np.pi/180)    
    
  return  np.real(spec), np.imag(spec), anguloOptimo


#%%

# PARA CORRERLO SOLO

# import matplotlib.pyplot as plt
# precision=1
# spec = real + 1j*imag  
# plt.figure(15433)
# plt.plot(ppmAxis, real)

# spec = spec[ppmAxis.argsort()] # ordeno basandome en ppmAxis
# ppmAxis = ppmAxis[ppmAxis.argsort()] # ordeno basandome en ppmAxis


# angle=np.arange(-180,180,precision)

# IntImagSpec = []
# IntRealSpec = []
# for i in range(angle.size):
#     Sp_try= spec*np.exp(-1j*angle[i]*np.pi/180)    
#     IntImagSpec.append(simps(np.imag(Sp_try),x=ppmAxis))
#     IntRealSpec.append(simps(np.real(Sp_try),x=ppmAxis))
    
# IntImagSpec = np.array(IntImagSpec)
# IntRealSpec = np.array(IntRealSpec)

# dataArray = np.array([np.abs(IntImagSpec), IntRealSpec, angle]).T
# dataArray = dataArray[dataArray[:, 0].argsort()] # ordeno basandome en la primera columna

# for i in range(IntImagSpec.size):
#   integralReal = dataArray[i,1]
#   if integralReal<0:
#     print("NO  ", dataArray[i,:])
#     continue
#   else:
#     print("SI  ", dataArray[i,:])
#     anguloOptimo = dataArray[i,2]
#     break
  
# spec = spec*np.exp(-1j*anguloOptimo*np.pi/180)    
# plt.figure(15433)
# plt.plot(ppmAxis, np.real(spec))

# plt.figure(1)
# plt.subplot(2,1,1)
# plt.plot(angle, IntImagSpec, label="imag")
# plt.legend()
# plt.subplot(2,1,2)
# plt.plot(angle, IntRealSpec, label="real")
# plt.legend()