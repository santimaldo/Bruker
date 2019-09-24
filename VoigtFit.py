# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:06:23 2019

@author: santi
"""
import numpy as np
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
        ===================COSAS A RESOLVER 24/09/2019==================
        ==    Implementar:                                            ==
        ==        Dar parametros iniciales                            ==
        ==        Fijar parametros o rangos                           ==
        ================================================================
    
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
    def __init__(self, x, y, Npicos=1, ajustar=True, **parametros_iniciales):
        
        self.x = x
        self.y = y
        self.Npicos = Npicos                        
        self.parametros_iniciales = parametros_iniciales
        self.modelo = None
        self.params = None
        self.ajuste = None
        
        # ajusto:
        self.generar_modelo()
        if ajustar:        
            self.fit()
    
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
            # chequeo si le di parametros iniciales, sino, usa los defaults
            if bool(self.parametros_iniciales):
                model_params = model.make_params(**self.parametros_iniciales)
            else:
                model_params = model.make_params(**default_params)
            
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
        self.params = params
    #--------------------------------------------------------------------------
        
    def fit(self):
        model = self.modelo            
        self.ajuste = model.fit(self.y, self.params, x=self.x)
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
    def plot(self):
        self.ajuste.plot()
