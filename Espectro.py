# coding=UTF-8
import nmrglue as ng
import numpy as np

class Espectro(object):

  """
  :version: 0
  :author: Santiago Agustin Maldonado Ochoa
  """

  """
    Clase Acqus. Contine informacion de los parametros usados en la medicion.


    Parameteros
    -----------
    path : str
        Es la carpeta donde se encuentran los datos de Bruker.
        Los atributos son los parametros mas usados para el analisis de datos.

    Atributos
    ----------
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

  """
  def __init__(self, path):
      size, real = ng.fileio.bruker.read_pdata_binary(path+'pdata/1/1r', big=False)

      self.size = int(size['FILE_SIZE'])
      self.real = real
      self.imag = ng.fileio.bruker.read_pdata_binary(path+'pdata/1/1i', big=False)[1]

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
