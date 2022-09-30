# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:47:43 2020

@author: santi
"""
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Exsy import *

from scipy.signal import savgol_filter
#%%


class Parametros(object):
    """
    defino esta clase para guardar los parametros que se ingresan a la funcion gauss
    """
    def __init__(self, centro):
        self.centro = centro
        self.Npeaks = 0
        self.A = []
        self.w = []
        self.Atot = 0
        
    def agregar(self, amp, w):
        self.A.append(amp)
        self.w.append(w)
        self.Npeaks += 1
        self.Atot = sum(self.A)

class Superficie(object):
    """
    aca voy guardando la superficie
    """
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y
        self.superficie = np.zeros_like(X)
        # este es el ultimo pico que se agrego
        self.ultimo = np.zeros_like(X)
        
    def agregar_pico(self, pico):
        """
        pico debe ser un array con la misma forma que X e Y
        """
        self.superficie = self.superficie + pico
        self.ultimo = pico
    
    def quitar_ultimo(self):
        """
        para quitar el ultimo que se agrego
        """
        self.superfice = self.superficie - ultimo
        self.ultimo = self.ultimo*0
        
def gauss2d(X, Y, par_x, par_y, angulo=-45):
    """
    par_x y par_y son objetos de tipo parametros
    """
    pi = np.pi       
    ang = pi/180*angulo
    gauss = np.zeros_like(X)
    
    xc = par_x.centro
    yc = par_y.centro
    
    # la amplitud es el promedio entre las ampl calculadas en cada dir.
    A = np.mean([par_x.Atot, par_x.Atot])
    
    for n in range(par_y.Npeaks):    
        
        # porcentaje del total:
        c_n = par_y.A[n]/par_y.Atot
        
        A_n = c_n * A
        
        wx2 = par_x.w[0]**2
        wy2 = par_y.w[n]**2
        cos2 = np.cos(ang)**2
        sen2 = np.sin(ang)**2
        
        # parametros para rotar la elipse
        # https://en.wikipedia.org/wiki/Gaussian_function
        a = cos2 /2/wx2 + sen2/2/wy2
        b = -np.sin(2*ang)/4/wx2 + np.sin(2*ang)/4/wy2
        c = sen2/2/wx2 + cos2/2/wy2        
        
        arg = a*(X-xc)**2 + 2*b*(X-xc)*(Y-yc) + c*(Y-yc)**2
        gauss = gauss + A_n * np.exp(-arg)  
    return gauss


#%%

path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY/25/"
savepath = "S:/Doctorado/Carbones/analisis/2020-01_CM7_EXSY/"
archivo_out = '1ms'

datos = DatosProcesados2D(path)

datos.espectro.ppmSelect2D([-8, 8])
# anulo todos los valores menores a 0:
datos.espectro.set_mask(1)


X = datos.espectro.ppmGrid_Dir
Y = datos.espectro.ppmGrid_Ind

#%%

superficie = Superficie(X,Y)


centrosx = [-3.0      , 1.9           , 2.87         , 4.51         ]
centrosy = [-3.0      , 1.9           , 2.87         , 4.48         ]
Ax =   [14112        , 6499          , 6000         , 81937        ]
wx =   [1.04384      , 0.44172       , 2.10966      , 0.3522       ]
Ay =   [(10543,3166) , (4141,8311)   , (4044,2004)  , (70894,9996)]
wy =   [(0.188,0.537), (0.696, 0.254), (0.340,0.244), (0.167,0.474)]


for n in range(len(centros)):
    par_x = Parametros(centrosx[n])
    par_y = Parametros(centrosy[n])
    
    par_x.agregar(Ax[n], wx[n])
    par_y.agregar(Ay[n][0], wy[n][0])
    par_y.agregar(Ay[n][1], wy[n][1])
    
    superficie.agregar_pico(gauss2d(X,Y, par_x, par_y))


surf = superficie.superficie

#%%
resta = datos.espectro.real - surf    
np.savetxt(savepath+'a'+'_diagonal.dat', surf)
np.savetxt(savepath+'a'+archivo_out+'-diagonal.dat', resta)
np.savetxt(savepath+'a_'+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
np.savetxt(savepath+'a_'+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)
    
    
N = 50
    
plt.figure(1)
plt.contour(X,Y,surf,N)
ax = plt.gca()
ax.invert_xaxis()
ax.invert_yaxis()

plt.figure(2)
plt.contour(X,Y,datos.espectro.real,N)
ax = plt.gca()
ax.invert_xaxis()
ax.invert_yaxis()

plt.figure(3)
plt.contour(X,Y,resta,N)
ax = plt.gca()
ax.invert_xaxis()
ax.invert_yaxis()


