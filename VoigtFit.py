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

        if fijar == []:
            self.anyfixed=False
            # si no se fijan parametros, se dejan todos los parametros como
            # variables. Esto es, se ajusta todo.
            for param in params:
                params[param].vary = True
        else:
            self.anyfixed=True

        
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
        
        #print(params['m1_center'].value, params['m1_center'].value)
        #print(params['m1_sigma'].value, params['m1_gamma'].value)
        self.ajuste = model.fit(self.y, params, x=self.x,
                                fitkws={'tol':1e-15})
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
    def plot_ajuste(self, reverse_xaxis=True):
        fig = self.ajuste.plot()
        ax = fig.gca()
        ii=0
        for component in self.componentes(self.x)[1]:
            ii+=1
            ax.plot(self.x, component, label=f'model {ii}')
        if reverse_xaxis:
            ax.set_xlim(np.max(self.x), np.min(self.x))
        return fig
    #--------------------------------------------------------------------------
    def plot_modelo(self):
        x = self.x
        y = self.modelo.eval(x=x, params=self.params)
        
        plt.figure(0)
        plt.plot(x,y)
        plt.show()
        

# -*- coding: utf-8 -*-
"""
Clase para ajuste de picos con PseudoVoigt usando lmfit.
Basado en la versión VoigtFit original.

@author: santi
"""


class PseudoVoigtFit:
    """
    Ajuste de espectros usando una suma de funciones PseudoVoigt.

    Attributes
    ----------
    x : array_like
        Eje x
    y : array_like
        Eje y
    Npicos : int
        Número de picos a ajustar
    modelo : lmfit.Model
        Modelo compuesto
    params : lmfit.Parameters
        Parámetros del modelo
    ajuste : lmfit.ModelResult
        Resultado del ajuste

    Methods
    -------
    generar_modelo():
        Genera el modelo compuesto con Npicos.
    generar_params():
        Asigna parámetros iniciales desde self.parametros_iniciales.
    fit(fijar):
        Realiza el ajuste, fijando parámetros si se pide.
    componentes(x):
        Devuelve (total, componentes) del ajuste evaluado en x.
    plot_ajuste():
        Grafica ajuste y componentes.
    plot_modelo():
        Grafica solo el modelo con los parámetros actuales.
    """

    def __init__(self, x, y, params=None, Npicos=1, ajustar=True, fijar=None, bounds=None, **parametros_iniciales):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.Npicos = Npicos
        self.parametros_iniciales = parametros_iniciales
        self.modelo = None
        self.params = params
        self.bounds = bounds
        self.ajuste = None

        if fijar is None:
            fijar = []

        # Si no se da un objeto params, lo inicializo
        NotInitParam = self.params is None

        self.generar_modelo()

        if bool(self.parametros_iniciales) and NotInitParam:
            self.generar_params()

        if ajustar:
            self.fit(fijar)
 
    # ------------------------------------------------------------------
    def generar_modelo(self):
        modelo_compuesto = None
        params = None
        x_min, x_max = np.min(self.x), np.max(self.x)
        x_range = x_max - x_min
        y_max = np.max(self.y)

        for i in range(self.Npicos):
            prefix_i = f"m{i+1}_"
            model = lm.models.PseudoVoigtModel(prefix=prefix_i)


            # Hints de parámetros
            model.set_param_hint("sigma", min=1e-3, max=x_range)
            model.set_param_hint("center", min=x_min, max=x_max)
            model.set_param_hint("height", min=1e-6, max=1.1 * y_max)
            model.set_param_hint("amplitude", min=1e-6)
            model.set_param_hint("fraction", min=0.0, max=1.0)

            
            #if i > 0:
            #    print("WARNING: sigma_i = sigma_1 for all peaks")
            #    model.set_param_hint('sigma', expr='m1_sigma')
            # Defaults aleatorios
            default_params = {
                prefix_i + "center": x_min + x_range * np.random.random(),
                prefix_i + "height": y_max * np.random.random(),
                prefix_i + "sigma": x_range / (2 * self.Npicos) * np.random.random(),
                prefix_i + "fraction": np.random.random(),
            }

            model_params = model.make_params(**default_params)

            # --- Update bounds si existen ---
            if self.bounds is not None:
                for pname, limits in self.bounds.items():
                    if "_" not in pname:
                        if pname in model_params:
                            model_params[pname].min = min(limits)
                            model_params[pname].max = max(limits)
                    else:
                        if pname in model_params:
                            model_params[pname].min = min(limits)
                            model_params[pname].max = max(limits)

            if self.params is None:
                if params is None:
                    params = model_params
                else:
                    params.update(model_params)

            if modelo_compuesto is None:
                modelo_compuesto = model
            else:
                modelo_compuesto += model

        self.modelo = modelo_compuesto
        if self.params is None:
            self.params = params

    # ------------------------------------------------------------------
    def generar_params(self):
        """
        Ejemplo:
            pf = PseudoVoigtFit(x, y, Npicos=2, center=[100,120], sigma=[5,10])
        """
        p_ini = self.parametros_iniciales
        for parametro, valores in p_ini.items():
            if not isinstance(valores, (list, tuple)):
                valores = [valores]
            for i, val in enumerate(valores):
                self.params[f"m{i+1}_{parametro}"].set(value=val)

    # ------------------------------------------------------------------
    def fit(self, fijar):
        params = self.params
        p_ini = self.parametros_iniciales

        if not isinstance(fijar, list):
            fijar = [fijar]

        if fijar == []:
            for param in params:
                params[param].vary = True

        else:
            for param_fijo in fijar:
                if "_" in param_fijo:
                    params[param_fijo].vary = False
                elif param_fijo in p_ini:
                    for i in range(len(p_ini[param_fijo])):
                        params[f"m{i+1}_{param_fijo}"].vary = False
                else:
                    for i in range(self.Npicos):
                        params[f"m{i+1}_{param_fijo}"].vary = False

        self.ajuste = self.modelo.fit(self.y, params, x=self.x)
        self.params = self.ajuste.params

    # ------------------------------------------------------------------
    def componentes(self, x):
        total = self.modelo.eval(x=x, params=self.params)
        comps = self.ajuste.eval_components(x=x)
        componentes = [comps[k] for k in sorted(comps.keys())]
        return total, componentes

    # ------------------------------------------------------------------
    def plot_ajuste(self, reverse_xaxis=True,
                    xlabel="x", ylabel="y"):
        fig = self.ajuste.plot()
        colors = ['r', 'b', 'g', 'orange', 'darkviolet', 'cyan', 'magenta', 'yellow', 'brown', 'pink']
        total = 0
        for ii, comp in enumerate(self.componentes(self.x)[1], start=1):
            fig.gca().plot(self.x, comp,color=colors[ii-1],label=f"model {ii}")
            total += comp
        fig.gca().plot(self.x, total, 'k', label='total fit')
        fig.gca().legend()
        if reverse_xaxis:
            fig.gca().set_xlim(np.max(self.x), np.min(self.x))
        fig.gca().set_xlabel(xlabel)
        fig.gca().set_ylabel(ylabel)
        return fig

    # ------------------------------------------------------------------
    def plot_modelo(self):
        y = self.modelo.eval(x=self.x, params=self.params)
        plt.figure()
        plt.plot(self.x, y, label="modelo")
        plt.legend()
        plt.show()
