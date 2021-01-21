# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:06:23 2019

@author: santi
"""
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm


class VoigtFit:
    """
    Class VoigtFit:
        Esta clase sirva para hacer ajuste de picos de los espectros, usando
        la funcion Voigt. Puedo dar la opcion de no ajustar, para solo generar
        el modelo.

    
    Attributes
    ----------
    
    x : Array_like
        Eje x
    
    y : Array_like
        Eje y
    
    Npicos : int
        Numero de picos a ajustar
        
    modelo : objeto de clase lmfit.Model
        Aqui se guarda el modelo que se ajusta
    
    params : objeto de clase lmfit.Parameters
        Guarda los parametros del modelo
    
    ajuste : objeto de clase lmfit.ModelResult
        Resultado del ajuste usando model.fit
        
    Methods
    -------
    generar_modelo():
        Como su nombre indica, genera el objeto Model. Lo que hace es agregar la
        cantidad de picos necesarios.
    
    fit():
        Realiza el ajuste utilizando Model.fit(). El resultado se guarda en el
        atributo ajuste, y el atributo params se modifica con los parametros
        correspondientes al mejor ajuste.

    componentes(x)
        Devuelve los datos ajustados. Hay que darle el eje x.
        return : total, componentes
            total : vector con la funcion de ajuste.
            componentes : lista con cada componente del ajuste, es decir, cada
                pico.
    
    plot():
        Es basicamente el metodo plot de la clase ModelResult. 
        Grafica el resultado y los residuos.
    
    """
    def __init__(self, x, y, params = None, Npicos=1, ajustar=True, fijar=[], **parametros_iniciales):
        
        self.x = x
        self.y = y
        self.Npicos = Npicos                        
        self.parametros_iniciales = parametros_iniciales
        self.modelo = None
        self.params = params
        self.ajuste = None
        
        if self.params is None:
            # esto es verdad si no le doy de entrada un OBJETO params
            NotInitParam = True
        else:
            NotInitParam = False
        
                                
        self.generar_modelo()
        # chequeo si le di parametros iniciales, en tal caso, se los aplico al
        # modelo. Notar que los parámetros que no fueron dados, seran los
        # aleatorios
        if bool(self.parametros_iniciales) and NotInitParam:
            p_ini = self.parametros_iniciales
            for parametro in p_ini:
                # si el parametro inicial es una lista, queda.
                if isinstance(p_ini[parametro], list):
                    continue
                # si tiene un solo elemento, lo transformo en lista
                elif type(p_ini[parametro]) in [int, float, str]:                    
                    p_ini[parametro] = [p_ini[parametro]]                
            self.generar_params()
        # ajusto:
        if ajustar:        
            self.fit(fijar)
    
    #----------------------------------Methods
    def generar_modelo(self):
        
        modelo_compuesto = None
        params = None
        x = self.x
        y = self.y
        Npicos = self.Npicos
        x_min = np.min(x)
        x_max = np.max(x)
        x_range = x_max - x_min
        y_max = np.max(y)
        
        for i in range(Npicos):
            prefix_i = f'm{i+1}_'
            model =  lm.models.VoigtModel(prefix=prefix_i)
            model.set_param_hint('sigma', min=1e-6, max=x_range)
            model.set_param_hint('gamma', min=1e-6, max=x_range)
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('height', min=1e-6, max=1.1*y_max)
            model.set_param_hint('amplitude', min=1e-6)            
            
            # default guess is horrible!! do not use guess()
            default_params = {
                prefix_i+'center': x_min + x_range/2 + x_range/2 * (np.random.random()-0.5),
                prefix_i+'height': y_max * np.random.random(),
                prefix_i+'sigma': x_range/2/Npicos * np.random.random(),
                prefix_i+'gamma': x_range/2/Npicos * np.random.random()
            }
            model_params = model.make_params(**default_params)
            
            if self.params is None:
                if params is None:
                    params = model_params
                else:
                    params.update(model_params)
                
            if modelo_compuesto is None:
                modelo_compuesto = model
            else:
                modelo_compuesto = modelo_compuesto + model                

        # return:                
        self.modelo = modelo_compuesto
        if self.params is None:
            self.params = params
    #--------------------------------------------------------------------------
    def generar_params(self):
        """
        La forma de darle los parametros inciales es, por ejemplo:
            vf = VoigtFit(Npicos=3, sigma=[10,20,20], center={246,260,270})
        """
        p_ini = self.parametros_iniciales        
        for parametro in p_ini:
            valores = p_ini[parametro]
            for i in range(len(valores)):
                # ojo, esto elimina los LIMITES dados en generar_modelo()
                self.params[f'm{i+1}_{parametro}'].set(value=valores[i])
        
        
        
    #--------------------------------------------------------------------------        
    def fit(self, fijar):
        model = self.modelo
        params = self.params
        p_ini = self.parametros_iniciales
        
        #si fijar no es una lista, lo convierto en una
        if not isinstance(fijar, list):
            fijar = [fijar]
        
        for param_fijo in fijar:
            # si hay un guion bajo, es porque solo se quiere fijar un pico.
            # ej: m2_center
            if '_' in param_fijo:
                params[f'{param_fijo}'].vary = False
            # si no se especifica en que pico, se fija ese parametro en todos 
            # los picos a los cuales se les da un valor inicial. Si ningun pico
            # tiene valor inicial, se fijan todos
            elif param_fijo in p_ini:
                for i in range(len(p_ini[param_fijo])):
                    params[f'm{i+1}_{param_fijo}'].vary = False
            elif param_fijo not in p_ini:
                for i in range(self.Npicos):
                    params[f'm{i+1}_{param_fijo}'].vary = False
        
        print(params['m1_center'].value, params['m1_center'].value)
        print(params['m1_sigma'].value, params['m1_gamma'].value)
        self.ajuste = model.fit(self.y, params, x=self.x)
        self.params = self.ajuste.params
                   
    #--------------------------------------------------------------------------    
    def componentes(self, x):
        total = self.modelo.eval(x=x, params=self.params)
        comps = self.ajuste.eval_components(x=x)
        componentes = []
        for i in range(len(comps)):
            componentes.append(comps[f'm{i+1}_'])

        return total, componentes
    #--------------------------------------------------------------------------
    def plot_ajuste(self):
        self.ajuste.plot()
        
    #--------------------------------------------------------------------------
    def plot_modelo(self):
        x = self.x
        y = self.modelo.eval(x=x, params=self.params)
        
        plt.figure(0)
        plt.plot(x,y)
        plt.show()
        

