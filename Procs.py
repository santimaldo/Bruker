# coding=UTF-8
import nmrglue as ng

class Procs(object):

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
      self.dic = ng.fileio.bruker.read_procs_file(path)['procs']
      self.ppm = self.dic["SF"]
      self.offset = self.dic["OFFSET"]
      self.FTsize = self.dic["FTSIZE"]
    
  def UnMetodo(self):
    """
    @return  :
    @author
    """
    pass
