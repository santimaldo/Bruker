# coding=UTF-8
import nmrglue as ng

class Acqus(object):

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
      self.dic = ng.bruker.read_acqus_file(path)['acqus']
      self.NS = int(self.dic["NS"])      
      self.PL1 = self.dic["PL"][1]
      self.P1 = self.dic["P"][1]
      self.RG = int(self.dic["RG"])
      self.SW = self.dic["SW"]

  def UnMetodo(self):
    """
    @return  :
    @author
    """
    pass
