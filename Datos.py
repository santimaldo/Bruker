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
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.

    Atributos
    ----------
    directorio : str
        En este atributo guardamos el nombre de la carpeta
    acqus : Acqus
        Es un objeto de la clase Acqus. Contiene los parametros con los cuales se realizo la medicion.
    procs : Procs
        Es un objeto de la clase Procs. Contiene los parametros con los cuales se procesaron los datos en TopSpin.
    pulseprog : PulseProg
        Es un objeto de la clase Procs. Contine informacion sobre el programa de pulso utilizado
    """
    def __init__(self, directorio):
        
        self.directorio = directorio
        self.acqus = Acqus(directorio)
        self.procs = Procs(directorio)
        self.pulseprog = PulseProg(directorio)

    def UnMetodo(self):
        """
        @return  :
        @author
        """
        pass

#----------------------------HERENCIAS------------------------------------------

class DatosProcesados(Datos):
    
    def __init__(self, directorio):
        Datos.__init__(self, directorio)
        self.espectro = Espectro()

    def set_espectro(self):
        """
        @return
        @author
        """
        path = self.directorio + 'pdata/1/1r'
        size, real = ng.fileio.bruker.read_pdata_binary(path, big=False)
        path = self.directorio + 'pdata/1/1i'
        null, imag = ng.fileio.bruker.read_pdata_binary(path, big=False)                
        
        self.espectro = Espectro()
        self.espectro.set_real(real)
        self.espectro.set_imag(imag)
        self.espectro.set_size(int(size['FILE_SIZE']))