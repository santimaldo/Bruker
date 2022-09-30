# coding=UTF-8
import nmrglue as ng

class Procs(object):

  """
  :version: 0
  :author: Santiago Agustin Maldonado Ochoa
  """

  """
    Clase Procs. Contine informacion de los parametros usados en el 
    pocesamiento de datos en TopSpin.
    Los atributos son los parametros mas usados para el analisis de datos.

    Parameteros
    -----------
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.
    p_dir : int
        Es la carpeta donde se encuentran los datos procesados.      

    Atributos
    ----------
    dic : dictionary
        Diccionario que contiene todos los parametros de procs
    FTsize : tuple
        Dimensiones del espectro
    ppm : float
        Numero de Hertz que representan 1 ppm, i.e frecuencia de resonancia
        dividido por 1e6.
    offset : float
        primer valor en el eje ppm.   
  """
  def __init__(self, directorio, p_dir, dim2 = False):
      path = directorio + 'pdata/' + str(p_dir) + '/'      
      
      procs = 'procs'
      if dim2:
          procs = 'proc2s'
          
      self.dic = ng.fileio.bruker.read_procs_file(path)[procs]
      self.FTsize = self.dic["FTSIZE"]
      self.ppm = self.dic["SF"]
      self.offset = self.dic["OFFSET"]
      self.phase = self.dic["PHC0"]
    
  def UnMetodo(self):
    """
    @return  :
    @author
    """
    pass
