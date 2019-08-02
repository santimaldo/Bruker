# coding=UTF-8

# librerias para la clase Datos
from Acqus import *
from Procs import *
from PulseProg import *
# librerias para la clase DatosProcesados
from Espectro import *

class Datos(object):

  """
  :version: 0
  :author: Santiago Agustin Maldonado Ochoa
  """

  """
    Clase Datos.

    Parameteros
    -----------
    path : str
        Es la carpeta donde se encuentran los datos de Bruker.

    Atributos
    ----------
    path : str
        En este atributo guardamos el nombre de la carpeta
    acqus : Acqus
        Es un objeto de la clase Acqus. Contiene los parametros con los cuales se realizo la medicion.
    procs : Procs
        Es un objeto de la clase Procs. Contiene los parametros con los cuales se procesaron los datos en TopSpin.
    pulseprog : PulseProg
        Es un objeto de la clase Procs. Contine informacion sobre el programa de pulso utilizado
  """
  def __init__(self, path):
    self.path = path
    self.acqus = Acqus(path)
    self.procs = Procs(path)
    self.pulseprog = PulseProg(path)

  def UnMetodo(self):
    """
    @return  :
    @author
    """
    pass

#----------------------------HERENCIAS------------------------------------------

class DatosProcesados(Datos):

  def __init__(self, path):
    Datos.__init__(self, path)
    self.espectro = Espectro(path)

  def UnMetodo(self):
    """
    @return  :
    @author
    """
    print('soy un metodo de DatosProcesados!')
