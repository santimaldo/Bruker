#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:33:48 2020

@author: santi
"""

import numpy as np
import matplotlib.pyplot as plt
from Espectro import *

#______________________________________________________________________________         
class Exsy(object):
    """
    Defino una clase para hacer el analisis de datos
    """
    
    def __init__(self, espectro, Npicos):
        """
        espectro es un objeto de la clase ESPECTRO.
        Npicos es el numero de picos diagonales
        """
        self.Npicos = Npicos
        self.spec = espectro.real
        self.ppmGridDir = espectro.ppmGrid_Dir
        self.ppmGridInd = espectro.ppmGrid_Ind     
        self.regiones = self.crear_regiones(Npicos)
        self.integrales = None
        
#------------------------------------------------------------------------------        
    def crear_regiones(self, Npicos):
        """
        crea la lista de listas vacias, donde guardaré los objetos REGION
        """
        regiones = []
        for n in range(Npicos):
            regiones.append([])
        for lista in regiones:
            for n in range(Npicos):
                lista.append(None)                
        return regiones
    
#------------------------------------------------------------------------------    
    def establecer_region(self, region, guess, delta_ppm=(0.6,0.6)):
        """
        este metodo establece la region de integracion para determinado pico.
        Inputs:
            region: tuple (n,m). Esto genera: self.regiones[n][m]
            guess: tuple(x,y). en que ppm se encuentra aproximadamente el pico
            delta_ppm: timple(delta_x, delta_y). cuanto se corre del guess para
                buscar el maximo
        """        
        # obtengo los indices del centro del pico.
        #xc, yc = self.encontrar_picos(guess, delta_ppm)
        xc,yc = guess
        x_index = find_nearest(self.ppmGridDir[0,:], xc)
        y_index = find_nearest(self.ppmGridInd[:,0], yc)
        # obtengo las coordenadas que determinan el rectangulo donde voy a
        # integrar.        
        x_lims, y_lims = self.establecer_limites(x_index, y_index)
        
        xi,xf = x_lims
        yi,yf = y_lims
        spec = self.spec[yi:yf, xi:xf]
        ppmGridDir = self.ppmGridDir[yi:yf, xi:xf]
        ppmGridInd = self.ppmGridInd[yi:yf, xi:xf]
        
        # dado que la region se definio con las derivadas, el pico está a lo sumo en el borde.
        # entonces en esta region busco el maximo y eso me define xc, yc.
        # estos son indices en la matriz reducida
        x_index = np.where(spec==np.max(spec))[1][0]
        y_index = np.where(spec==np.max(spec))[0][0]
        
        # ahora busco cuanto vale la coordenadas en esos indices para luego encontrar
        # los indices en la matriz completa
        xc = ppmGridDir[0,x_index]
        yc = ppmGridInd[y_index,0]      
        # estos son indices en la matriz completa
        x_index = find_nearest(self.ppmGridDir[0,:], xc)
        y_index = find_nearest(self.ppmGridInd[:,0], yc)        
                            
        # ahora repito el procedimiento para determinar la region, a partir de un nuevo xc, yc
        x_lims, y_lims = self.establecer_limites(x_index, y_index)
        
        print(x_lims, y_lims)
        
        xi,xf = x_lims
        yi,yf = y_lims
        spec = self.spec[yi:yf, xi:xf]
        ppmGridDir = self.ppmGridDir[yi:yf, xi:xf]
        ppmGridInd = self.ppmGridInd[yi:yf, xi:xf]
                
        n, m = region
        self.regiones[n][m] = Region(ppmGridDir, ppmGridInd, spec)
        return 0
        
#------------------------------------------------------------------------------    
    def encontrar_picos(self, guess, delta_ppm):
        """
        Esta funcion busca el maximo de spec alrededeor de un punto x,y, en un
        cuadrado de 2*delta_x por 2*delta_y
        
        Inputs:
            region: tuple (n,m). Esto genera: self.regiones[n][m]
            guess: tuple(x,y). en que ppm se encuentra aproximadamente el pico            
            
            #delta_x, delta_y sirven para definir nx, ny, esto es: cuantos pixels 
            #definir el lugar para buscar el maximo
            delta_x : ppm hacia izquierda y derecha
            delta_y : ppm hacia arriba y abajo
        """
        # La idea es esta:
        # dado el delta en cada direccion, defino una region cuadrada donde
        # busco el maximo. Luego, con los indices del maximo encuentro el pico.
        # es importante que el guess sea bueno! 
        ppmDir = self.ppmGridDir[0,:]
        ppmInd = self.ppmGridInd[:,0]
        
        x,y = guess        
        x_index = find_nearest(ppmDir, x)
        y_index = find_nearest(ppmInd, y)
        
        # cuantos ppm son un paso en cada direccion
        stepDir = np.abs(ppmDir[1]-ppmDir[0])
        stepInd = np.abs(ppmInd[1]-ppmInd[0])            
        nx = int(delta_ppm[0]/stepDir)
        ny = int(delta_ppm[1]/stepInd)
        
        spec_reduced = self.spec[y_index-ny:y_index+ny, x_index-nx:x_index+nx]
        ppmDir_reduced = ppmDir[x_index-nx:x_index+nx]
        ppmInd_reduced = ppmInd[y_index-ny:y_index+ny]
        
        maximo = spec_reduced.max()
        y, x = np.where(spec_reduced==maximo)
        
        x = ppmDir_reduced[x[0]]
        y = ppmInd_reduced[y[0]]
        
        x_index = find_nearest(ppmDir, x)
        y_index = find_nearest(ppmInd, y)
        
        #plt.contourf(ppmDir_reduced,ppmInd_reduced,spec_reduced)
        return x_index, y_index
#------------------------------------------------------------------------------
    def establecer_limites(self, x_index, y_index):
        """
        esta funcion debe devolver dos tuplas (xi,xf), (yi,yf)
        """
#        x_index = int(x_index)
#        y_index = int(y_index)
        
        spec_x = self.spec[y_index, :]
        spec_y = self.spec[:, x_index]

        dspec_x = np.diff(spec_x)
        dspec_y = np.diff(spec_y)
        # la derivada tiene un pubnto menos, defino nuevos ejes
        dppmDir = self.ppmGridDir[0, range(dspec_x.size)]
        dppmInd = self.ppmGridInd[range(dspec_y.size), 0]
        
        x_lims = self.find_zeros(dppmDir, dspec_x, x_index)
        y_lims = self.find_zeros(dppmInd, dspec_y, y_index)
        
        
        # grafico para testear:
        plt.figure(4291)
        ppmDir = self.ppmGridDir[0,:]
        ppmInd = self.ppmGridInd[:,0]
        
        xi,xf = x_lims
        yi,yf = y_lims
        plt.subplot(2,1,1)
        plt.plot(ppmDir,spec_x, 'b--')
        plt.plot(ppmDir[xi:xf],spec_x[xi:xf], 'r')
        plt.subplot(2,1,2)
        plt.plot(ppmInd,spec_y, 'b--')
        plt.plot(ppmInd[yi:yf],spec_y[yi:yf], 'r')
        ax = plt.gca()
        ax.invert_xaxis()
        
        
        return x_lims, y_lims
#------------------------------------------------------------------------------
    def find_zeros(self, x_axis, array, inicio_index):        
        """
        el objetivo de esta funcion es darle la derivada del spec, y que desde el 
        punto inicio_index, recorra el array hasta que llegue a cero, o bien,
        que se estanque en algun valor constante
        """
        delta_ppm = 0.2
        step = np.abs(x_axis[1]-x_axis[0])
        delta_n = int(delta_ppm/2/step)
                
        mean_x = []
        mean_array = []
        std_array = []
        maximo = abs(array).max()
        
        index = inicio_index + 4*delta_n
        intervalo = array[index-delta_n : index+delta_n]
        new = np.mean(intervalo)
        condition = True
        while condition:            
            old = new
            sgn_old = np.sign(old)
            std = np.std(intervalo)
            # guardo los valores
            mean_x.append(x_axis[index])
            mean_array.append(old)
            std_array.append(std)
            # doy un paso
            index = int(index + 2*delta_n)
            if index + 2*delta_n < array.size:
                intervalo = array[index-delta_n : index+delta_n]    
            else:
                break
            new = np.mean(intervalo)
            sgn_new = np.sign(new)              
#            cond1 = sgn_new==sgn_old
#             # cond2: es falsa cuando la desviacion es chica y esta cerca de cero)
#            cond2 = not( abs(std)<0.01*maximo and abs(new)<0.01*maximo)
#            condition = cond1 and cond2
            condition = sgn_new==sgn_old
        
        index_fin = index - delta_n
            
        index = inicio_index - 4*delta_n
        intervalo = array[index-delta_n : index+delta_n]
        new = np.mean(intervalo)
        condition = True
        while condition:
            old = new
            sgn_old = np.sign(old)
            std = np.std(intervalo)
            # guardo los valores
            mean_x.append(x_axis[index])
            mean_array.append(old)
            std_array.append(std)
            # doy un paso
            index = int(index - 2*delta_n)
            if index - 2*delta_n > 0:
                intervalo = array[index-delta_n : index+delta_n]    
            else:
                break
            new = np.mean(intervalo)
            sgn_new = np.sign(new)
#            cond1 = sgn_new==sgn_old
#            cond2 = not( abs(std)<0.01*maximo and abs(new)<0.01*maximo)
#            condition = cond1 and cond2
            condition = sgn_new==sgn_old
            
        index_ini = index + delta_n
        
        
        
        # grafico para testear:
        plt.figure(432)
        plt.plot(x_axis,array*0, 'k--')
        plt.plot(x_axis,array)
        plt.plot(mean_x, mean_array, 'o')
        ax = plt.gca()
        ax.invert_xaxis()
        
        plt.figure(431)
        plt.plot(mean_x, np.zeros_like(mean_x))
        plt.plot(mean_x, std_array, 'o')
        ax = plt.gca()
        ax.invert_xaxis()
        
        
        return index_ini, index_fin
#------------------------------------------------------------------------------
    def graficar_regiones(self, ax = None, Ncontour=30):

        # primero el contour del espectro entero        
        # plt.figure(int(fignum))
        if ax is None:
              ax = plt.subplot(1,1,1)
              
        ax.contour(self.ppmGridDir, self.ppmGridInd, self.spec, Ncontour, cmap='binary')
        

        
        for lista in self.regiones:
            for region in lista:        
                if not region is None:
                    x = region.ppmGridDir
                    y = region.ppmGridInd
                    spec = region.spec
                    ax.contourf(x,y,spec, cmap='jet')
                    
        ax = plt.gca()
        ax.invert_yaxis()
        ax.invert_xaxis()
        return 0
#------------------------------------------------------------------------------
    def integrar_regiones(self):
        integrales = np.zeros([self.Npicos,self.Npicos])
        for i in range(self.Npicos):
            for j in range(self.Npicos):        
                try:
                    region = self.regiones[i][j]
                    integrales[i,j] = region.integrar()
                except:
                    continue
        self.integrales = integrales
        return integrales      
#____________________________________________________________________Fin_Exsy__
        
    
    
#______________________________________________________________________________    
class Region(object):
    """
    objeto con la region asociada a cada pico
    """
    def __init__(self, ppmGridDir, ppmGridInd, spec):
        self.spec = spec
        self.ppmGridDir = ppmGridDir
        self.ppmGridInd = ppmGridInd
        
        self.ppmDir = ppmGridDir[0,:]
        self.ppmInd = ppmGridInd[:,0]
        self.integral = None
    
    def plot_region(self):        
        plt.figure(1)
        plt.contourf(self.ppmGridDir, self.ppmGridInd, self.spec, cmap='jet')
        plt.show
        
    def integrar(self):
        """
        integracion de superficie
        """        
        x = self.ppmDir
        y = self.ppmInd
        I = np.trapz(self.spec, x=y, axis=0)
        I = np.trapz(I, x=x)
        self.integral = I
        return I
#______________________________________________________________________________


    
#____________________________FUNCIONES_________________________________________
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