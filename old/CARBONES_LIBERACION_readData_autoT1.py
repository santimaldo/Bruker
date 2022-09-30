# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 11:42:36 2018

@author: santi


LA CONCLUSION DE ESTO ES QUE EL T1 DEL AGUA SE MANTIENE CONSTANTE
"""
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from scipy.optimize import curve_fit
from scipy.integrate import simps


def ExpDec1(x, t1, A, y0):
	return A*np.exp(-x/t1)+y0

def ExpDec2(x, t1, t2, A1, A2, y0):
	return A1*np.exp(-x/t1)+A2*np.exp(-x/t2)+y0


#-------------------------- INPUTS --------------------------------------------
path0 = 'S:/Doctorado/Carbones/300MHz/2019-08-27_Carbones_Liberacion_CM7/'
#path0 = 'S:/Doctorado/Carbones/300MHz/2019-08-26_Carbones_Liberacion_CM12/'

path1 ='S:/Doctorado/Carbones/analisis/2019-08_Liberacion/CM7/T1/'
#------------------------------------------------------------------------------

#carpetas = [2,6,10,14,18,22,26,30,34,38,42,46]
#tiempos = [0, 0.5, 1, 2, 4, 6, 10, 15, 20, 25, 30, 280] # minutos


tiempos = [0,0.5,1,2,4,8,13,40,61,80]
# tiempos= [0]
carpetas = [(i*4)+2 for i in range(len(tiempos))] # agua



# carpetas = [38]
# tiempos = [80] # horas
unidad_de_tiempo = 'h' # para usar en el nombre de los archivos

monoexponencial = False
T1list = []

# frecuencias de integracion
ppmi = -20
ppmf =  20


for j in range(len(carpetas)):
    carpeta = carpetas[j]
    path = f"{path0}/{carpeta}/"
    archivo_out = 't_' + str(tiempos[j]) + unidad_de_tiempo
        
    #-------------------------- EXTRAER FACTOR VD----------------------------------
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
    vdlist = np.loadtxt(f"{path}/vdlist")          
    tauAxis = vdlist*factor_vd
    #------------------------------------------------------------------------------
    
    
    
    #-------------------------- IMPORTAR ------------------------------------------
    datos = DatosProcesados2D(f"{path0}/{carpeta}/")
    # datos.espectro.ppmSelect2D([-8, 8])
  
    spec = datos.espectro.real    
    ppmAxis = datos.espectro.ppmAxis
    
    T1data = np.zeros_like(spec[:,0])
    condicion = np.abs(ppmAxis)<20
    for jj in range(T1data.size):
      spec_tau = spec[jj,:]
      spec_tau = spec_tau[condicion]
      T1data[jj] = simps(spec_tau, x=-ppmAxis[condicion])
    T1data = T1data/np.max(np.abs(T1data))           
      
    tauAxis = tauAxis[0:T1data.size]    
    #------------------------------------------------------------------------------
    
    
    
    #-------------------------- AJUSTAR -------------------------------------------
    #guess = 185 * np.exp(j/20.0) -144.
    # guess = 200
    # popt, pcov = curve_fit(exponenial_func, x, y, p0=(1., guess, -1.),absolute_sigma=True)
    #guess = 185 * np.exp(j/20.0) -144.
    # t = np.linspace(tau[0], tau[-1], 100)    
    
    # if monoexponencial:
    #     try:
    #         popt, pcov = curve_fit(ExpDec1, tau, S,absolute_sigma=True)
    #         T1 = popt[0]; A = popt[1]; y0 = popt[2]
    #         ajuste = ExpDec1(t,T1,A,y0)
    #         print('T1= '+ str(T1))
    #         T1list.append(T1)
    #         no_error = True            
    #     except:
    #         print('Ocurrio un error en el ajuste monoexponencial')
    #         no_error = False
    #         T1list.append(0)
    # else:            
    # try:
    #   popt, pcov = curve_fit(ExpDec2, tauAxis, T1data, absolute_sigma=True)
    #   T1a = popt[0]; T1b = popt[1]; Aa = popt[2]; Ab = popt[3]; y0 = popt[4]
    #   ajuste = ExpDec2(tauAxis,T1a, T1b, Aa, Ab, y0)    
    #   no_error = True
    #   print('T1a= '+ str(T1a))
    #   print('T1b= '+ str(T1b))                        
    # except:
    #   no_error = False
    #------------------------------------------------------------------------------
    
    
    #-------------------------- GRAFICO -------------------------------------------
    plt.figure(j)
    plt.plot(tauAxis,T1data ,'o', label=tiempos[j])
    if no_error:
      plt.plot(tauAxis, ajuste, '-')
    
    # #------------------------------------------------------------------------------
    
    # #-------------------------- GUARDAR -------------------------------------------
    # # np.savetxt(path1+archivo_out+'.dat', T1data)
    
    # # creando el archivo .info
    # with open(path + 'pdata/1/title', 'rt') as f:
    #     titulo = f.read()
    
    # with open(path + 'pdata/1/outd', 'rt') as f:    
    #     for line in f:
    #       if line.lstrip().startswith('#'):
    #         continue
    #       elif line.lstrip().startswith('$'):    
    #         fecha = line.rstrip('\n').lstrip('$$')
    #         print(fecha)
    #         break
    
    
    # with open(path1+archivo_out+'.info', 'w') as f:
    
        
    #     f.write('datos:\n"'   +path+'"\n')
    #     f.write('factorvd:\n'+str(factor_vd)+'\n')
    #     f.write('fecha y hora de medicion:\n'   +fecha+'\n')
    #     f.write('\n')
    #     f.write('titulo:\n'+titulo)
  
    #------------------------------------------------------------------------------
    j += 1
    
plt.legend()    
plt.show    
