# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:06:23 2019

@author: santi
"""
import numpy as np
import lmfit as lm


class VoigtFit:
    """
    DOCUMENTACION
    """
    def __init__(self, x, y, Npicos=1,**parametros_iniciales):
        
        self.x = x
        self.y = y
        self.Npicos = Npicos                        
        self.parametros_iniciales = parametros_iniciales
        self.modelo = None
        self.params = None
        self.generar_modelo()

    def UnMetodo(self):
        """
        @return  :
        @author
        """
        pass
    
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
            # if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
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