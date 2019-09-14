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
        self.pulseprog = PulseProg(directorio)

    def UnMetodo(self):
        """
        @return  :
        @author
        """
        pass

#----------------------------HERENCIAS------------------------------------------        
class DatosProcesados(Datos):
    
    def __init__(self, directorio, pdata=1):
        Datos.__init__(self, directorio)        
        self.p_dir = p_dir
        self.procs = Procs(directorio, p_dir)
        self.espectro = Espectro()
        self.set_espectro()
        
        
        
    def set_espectro(self):
        """
        @return
        @author
        """            
        path = self.directorio + 'pdata/' + str(self.pdata) + '/'
        dic, data = ng.fileio.bruker.read_pdata(path)
        real = np.real(data)
        imag = np.imag(data)
                        
        self.espectro.set_real(real)
        self.espectro.set_imag(imag)
        self.espectro.set_size(real.shape)
        self.espectro.set_ppmAxis(self.crear_ppmAxis())

    def crear_ppmAxis(self):
        ppm = self.procs.ppm
        offset = self.procs.offset
        FTsize = self.procs.FTsize
        SpectralWidth = self.acqus.SW
        
        ultimo_ppm = offset - SpectralWidth
        ppmAxis = np.linspace(offset, ultimo_ppm, FTsize)
        return ppmAxis 

        