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
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.
        Los atributos son los parametros mas usados para el analisis de datos.

    Atributos
    ----------
    dic : dictionary
        Diccionario que contiene todos los parametros de acqus
    NS : int
        Numero de Scanes
    P1 : float
        Potencia del primer pulso (dB)
    PL1 : float
        Duracion del primer pulso (us).    
    RG : int
        Ganancia del Receptor        
    SW : float
        Ancho espectral en ppm.

  """
  def __init__(self, directorio):
      self.dic = ng.bruker.read_acqus_file(directorio)['acqus']
      self.NS = int(self.dic["NS"])      
      self.P1 = self.dic["P"][1]
      self.PL1 = self.dic["PL"][1]
      self.RG = int(self.dic["RG"])
      self.SW = self.dic["SW"]

  def UnMetodo(self):
    """
    @return  :
    @author
    """
    pass
