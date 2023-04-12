# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:09:33 2022

@author: Santi

Extrae multiples espectros adquiridos en Bruker Avance II
"""

# import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from Datos import *
from Espectro import autophase

#------------------------------------
path = "S:/PosDoc/Gicerol-agua/116MHz/2022-09-29_Diff_Silica_Agua-Glicerol-LiCl/"
savepath = 'S:/temp/'
expnums = np.arange(20,34)
muestra = "M10_Q3_50pc"
nucleo = "7Li"
#---------------------------------
path = "S:/PosDoc/Gicerol-agua/116MHz/2022-10-03_Diff_Silica_Agua-Glicerol-LiCl/"
savepath = 'S:/temp/'
expnums = np.arange(30,42)
expnums = np.arange(60,67)
muestra = "M10_Q3_50pc"
nucleo = "7Li"

ppmRange = [-2,0.5]
gplist = []
intensidades = []
intensidadesFID = []
fig, axs = plt.subplots(num=1,nrows=1, ncols=2 )  # create figure & 1 axis
for nn in range(len(expnums)):
  expn = expnums[nn]  
  print(expn)
  # extraigo:
  datos = DatosProcesados(f'{path}/{expn}/')
  gplist.append(datos.acqus.gp)
  datos.espectro.ppmSelect(ppmRange)
  re = datos.espectro.real
  im = datos.espectro.imag
  ppmAxis = datos.espectro.ppmAxis

  integral = scipy.integrate.simps(re, x=-ppmAxis)  
  intensidades.append(integral)
  
  ######### calculo FID
  datos.set_fid()
  timeAxis = datos.fid.timeAxis
  fid = datos.fid.real
  fid = np.abs(datos.fid.real + 1j*datos.fid.imag)
  npts=20
  intensidadFID = np.abs(np.sum(fid[0:npts]))
  intensidadesFID.append(intensidadFID)
  

  # guardo:
  # header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
  # dataexport = np.array([ppmAxis, re, im, re, im]).T
  # filename = f'{savepath}/{nucleo}_{muestra}.dat'
  # # np.savetxt(filename, dataexport, header=header)

  # grafico para ver:
  print('graficando...', nucleo, muestra)
  axs[0].plot(ppmAxis, re, linewidth=2)
  axs[0].set_xlabel(f"{nucleo} NMR Shift [ppm]")
  axs[0].set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
  axs[1].plot(timeAxis, fid, 'o-', linewidth=2)
  axs[1].set_xlabel(f"time [ms]")
  


  # # grafico para guardar:
  # fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
  # ax.plot(ppmAxis, re, linewidth=2)
  # ax.set_title(muestra)
  # ax.set_xlabel(f"{nucleo} NMR Shift [ppm]")
  # ax.set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
  # filename = f'{savepath}/{nucleo}_{muestra}.png'
  # fig.savefig(filename)   # save the figure to file
  # plt.close(fig)    # close the figure window

#%%
# calculo bvalue:
# gplist = np.loadtxt(f"{path}gp_list.dat")[:,1]
gplist = np.array(gplist)
bigDelta = 10e-3 # s
delta    = 5e-3 # s
gamma = 103.962e6  #  rad/(s*T) ---> 7Li
g0 = 12 # T/m

S = np.array(intensidades)
# S = np.array(intensidadesFID)


bvalue = (gamma*gplist/100*g0*delta)**2*(bigDelta-delta/3) * 1e-9
plt.figure(75311)
plt.semilogy(bvalue, S/S[0], 'o')



plt.show()
