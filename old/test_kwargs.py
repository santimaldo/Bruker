# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:27:59 2019

@author: santi
"""

class Clase(object):
    
    def __init__(self, **parametros):        
        self.parametros = parametros
        



a = Clase(arg1=1, arg2=2)
b = Clase()


print(a.parametros)
print(b.parametros)
        