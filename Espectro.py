# coding=UTF-8
import nmrglue as ng
import numpy as np

class Espectro(object):

    """
    :version: 0
    :author: Santiago Agustin Maldonado Ochoa
    """

    """
    Clase Espectro. Es la transformada de fourier de la fid. Estos objetos
    pueden provenir de un procesamiento en TopSpin o bien de procesar una fid
    externamente


    Parameteros
    -----------
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.
        
    Atributos
    ----------
    Los atributos son los parametros mas usados para el analisis de datos.

    dic : dictionary
        Diccionario que contiene todos los parametros de acqus
    NS : int
        Numero de Scanes
    RG : int
        Ganancia del Receptor
    PL1 : float
        Duracion del primer pulso (us).
    P1 : float
        Potencia del primer pulso (dB)

    Medotos
    -------
    
    abs():
          return: array
          devuelve el modulo del espectro.
          
    complex():
          return: array
          devuelve el expectro como numero complejo.

    ppmSelect(ppmRange):
          Redefine todos los parametros pero con el eje ppm recortado.
          
          parametros: list.
                ppmRange = [ppmInicial, ppmFinal]
    """
    def __init__(self):
        
        self.real = None
        self.imag = None
        self.ppmAxis = None
        self.size = None
        
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
  
    def ppmSelect(self):
        pass
