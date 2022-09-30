# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:00:55 2020

@author: santi
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter

from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion


#------------------------------------------------------------------------------
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
#------------------------------------------------------------------------------
def derivar(array, ws=151, po=3):
    """
    esta funcion deriva y suaviza
    ws: window size
    po: polynomial order
    """
    d = np.diff(array)
    d = savgol_filter(d, ws, po)
    
    # chequeo que no se muy grande el ws
    cond1 = np.abs(d[0])>=(np.abs(d).max()/20)
    cond2 = np.abs(d[-1])>=(np.abs(d).max()/20)
    while cond1 or cond2:               
        ws = ws-20
        d = np.diff(array)
        d = savgol_filter(d, ws, po)
        cond1 = np.abs(d[0])>=(np.abs(d).max()/10)
        cond2 = np.abs(d[-1])>=(np.abs(d).max()/10)
    return d
#------------------------------------------------------------------------------
def find_peak(ppmDir, ppmInd, spec, x, y, delta_x=2.5, delta_y=2.5):
    """
    Esta funcion busca el maximo de spec alrededeor de un punto x,y, en un
    cuadrado de 2*delta_x por 2*delta_y
    
    delta_x, delta_y sirven para definir nx, ny, esto es: cuantos pixels 
    definir el lugar para buscar el maximo
    delta_x : ppm hacia izquierda y derecha
    delta_y : ppm hacia arriba y abajo
    """
    x_index = find_nearest(ppmDir, x)
    y_index = find_nearest(ppmInd, y)
    
    # cuantos ppm son un paso en cada direccion
    stepDir = np.abs(ppmDir[1]-ppmDir[0])
    stepInd = np.abs(ppmInd[1]-ppmInd[0])
        
    nx = int(delta_x/stepDir)
    ny = int(delta_y/stepInd)
    spec_reduced = spec[y_index-ny:y_index+ny, x_index-nx:x_index+nx]
    ppmDir_reduced = ppmDir[x_index-nx:x_index+nx]
    ppmInd_reduced = ppmInd[y_index-ny:y_index+ny]
    
    maximo = spec_reduced.max()
    yy, xx = np.where(spec_reduced==maximo)
    
    x = ppmDir_reduced[xx[0]]
    y = ppmInd_reduced[yy[0]]
    
    x_index = find_nearest(ppmDir, x)
    y_index = find_nearest(ppmInd, y)
    
    plt.contourf(ppmDir_reduced,ppmInd_reduced,spec_reduced)
    
    return x_index, y_index
#------------------------------------------------------------------------------
def find_zeros(array, inicio_index, x_axis=None, N=5):
    """
    el objetivo de esta funcion es darle la derivada del spec, y que desde el 
    punto inicio_index, recorra el array hasta que llegue a cero
    """
    
    # cerca del centro pude pasar varias veces por el cero. A izquiuerda y 
    # derecha debe tener distintos signos. busco los puntos donde esto ocurre y
    # desde ahi arranco.
    sgn_left  = np.sign(array[inicio_index-N])
    sgn_right = np.sign(array[inicio_index+N])
    while sgn_left == sgn_right:
        sgn_left  = np.sign(array[inicio_index-N])
        sgn_right = np.sign(array[inicio_index+N])
        N += 1        
        
    n = inicio_index + N # no arranco desde el max para evitar errores de cambio de signo
    sgn_old = 1
    sgn_new = 1
    while sgn_old == sgn_new:
        old = array[n]
        sgn_old = np.sign(old)
        n += 1
        new = array[n]
        sgn_new = np.sign(new)
        #si no llega a cero pero es constante en variaciones de 1%, corta
        #if np.abs((new-old)/old)<0.01:
        #    break        
    fin_index = n    
    n = inicio_index - N # no arranco desde el max para evitar errores de cambio de signo
    sgn_old = 1
    sgn_new = 1
    while sgn_old == sgn_new and n>0:
        sgn_old = np.sign(array[n])
        n -= 1
        sgn_new = np.sign(array[n])
    ini_index=n
    
    if not x_axis is None:
        plt.figure(432)
        plt.plot(x_axis,array*0)
        plt.plot(x_axis,array)
        plt.plot(x_axis[ini_index:fin_index],array[ini_index:fin_index])
        plt.plot(x_axis[inicio_index+N],array[inicio_index+N], 'o')
        plt.plot(x_axis[inicio_index],array[inicio_index], 'o')
        plt.plot(x_axis[inicio_index-N],array[inicio_index-N], 'o')
    
    return ini_index, fin_index
#------------------------------------------------------------------------------
def integrar(x, y, surf):
    """
    integracion de superficie
    x, y array 1d de coordenadas.
    """        
    I = np.trapz(surf, x=y, axis=0)
    I = np.trapz(I, x=x)
    return I
#------------------------------------------------------------------------------

        
plt.figure(100)
plt.contour(ppm_x, ppm_y, spec, 25,  cmap='jet', vmax=5000000)
plt.show()


x = -3.0
y = 1.9
x_index, y_index = find_peak(ppmDir, ppmInd, spec, x, y, delta_x=0.8, delta_y=0.8)

spec_x = spec[y_index, :]
spec_y = spec[:, x_index]

#plt.figure(2)
#plt.plot(ppmDir, spec_x)
#plt.plot(ppmInd, spec_y)

dspec_x = derivar(spec_x)
dspec_y = derivar(spec_y)
dppmDir = ppmDir[range(dspec_x.size)]
dppmInd = ppmInd[range(dspec_y.size)]

#plt.figure(3)
#plt.plot(dppmDir, dspec_x)
#plt.plot(dppmInd, dspec_y)


ini_index_x, fin_index_x = find_zeros(dspec_x, x_index, x_axis=dppmDir)
plt.figure(4)
plt.plot(ppmDir, spec_x)
plt.plot(ppmDir[ini_index_x:fin_index_x], spec_x[ini_index_x:fin_index_x])


ini_index_y, fin_index_y = find_zeros(dspec_y, y_index, x_axis=dppmInd)
#plt.figure(5)
##plt.plot(ppmDir, spec_x)
#plt.plot(ppmInd[ini_index_y:fin_index_y], spec_y[ini_index_y:fin_index_y])
ini_index_y = 56

plt.figure(200)
plt.title('regiones a integrar')
plt.contour(ppm_x, ppm_y, spec, 25,  cmap='jet', vmax=5000000)
#slices:
slice_x = ppm_x[ini_index_y:fin_index_y,ini_index_x:fin_index_x]
slice_y = ppm_y[ini_index_y:fin_index_y,ini_index_x:fin_index_x]
slice_spec = spec[ini_index_y:fin_index_y,ini_index_x:fin_index_x]
plt.contourf(slice_x, slice_y, slice_spec,  cmap='jet')
plt.show()


slice_ppmDir = ppmDir[ini_index_x:fin_index_x]
slice_ppmInd = ppmInd[ini_index_y:fin_index_y]
I = integrar(slice_ppmDir, slice_ppmInd, slice_spec)
print(I)





xi_32 = ini_index_x
xf_32 = fin_index_x
yi_32 = ini_index_y
yf_32 = fin_index_y