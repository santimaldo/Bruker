# coding=UTF-8
import nmrglue as ng
import numpy as np

class Fid(object):

    """
    :version: 0
    :author: Santiago Agustin Maldonado Ochoa
    """

    """
    Clase Fid. Es el valor de la senal de RMN. Puede ser 1D o 2D.


    Parameteros
    -----------
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.
        Los atributos son los parametros mas usados para el analisis de datos.

    Atributos
    ----------
    atributo : tipo
        Explicacion de que es.

    """
    def __init__(self):
        
        self.real = 0
        self.imag = 0
        self.size = 1
        
    def set_real(self, real):
        self.real = real
    
    def set_imag(self, imag):
        self.imag = imag
        
    def set_size(self, size):
        self.size = size
    
    def set_ppmAxis(self, ppmAxis):
        self.ppmAxis = ppmAxis
        
    def abs(self):
        """
        @return  :
        @author
        """
        complex = self.real + 1j * self.imag
        return np.abs(complex)

    def complex(self):
        """
        @return  :
        @author
        """
        complex = self.real + 1j * self.imag
        return complex
