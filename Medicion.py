# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 18:53:52 2019

@author: santi
"""

class Medicion(object):
    """
    Clase Medicion
    --------------
   
    Es la clase de la cual se heredan las demas clases
   
    Atributos
    ---------
    directorio : string
        ubicacion del conjunto de archivos de la medicion.
        
    Metodos
    -------
    get_directorio()
        return : string
   
    """
    def __init__(self, directorio):
        self.directorio = directorio # string

    def get_directorio(self, directorio):
        return self.directorio