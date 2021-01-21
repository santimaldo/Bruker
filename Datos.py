# coding=UTF-8

# librerias para la clase Datos
from Acqus import *
from Procs import *
from PulseProg import *
# librerias para la clase DatosProcesados
from Espectro import *
from Fid import *

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

#----------------------------HERENCIAS-----------------------------------------
class DatosCrudos(Datos):
    def __init__(self, directorio):
        Datos.__init__(self, directorio)
        self.fid = Fid()
        self.set_fid()
    
    def set_fid(self):
        """
        @return
        @author
        """                    
        null, data = ng.fileio.bruker.read(self.directorio)                
        real = np.real(data)
        imag = np.imag(data)
        
        # elimino los primeros puntos del tiempo muerto      
        try:
          real = real[:,69:]
          imag = imag[:,69:]
        except:
          real = real[69:]
          imag = imag[69:]
        # creo el eje temporal
        timeAxis = np.arange(np.size(real)) * self.acqus.DW
        
                        
        self.fid.set_real(real)
        self.fid.set_imag(imag)
        self.fid.set_timeAxis(timeAxis)
        self.fid.set_size(real.shape)
        
        
class DatosProcesados(Datos):
    
    def __init__(self, directorio, p_dir=1):
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
        path = self.directorio + 'pdata/' + str(self.p_dir) + '/'                
        null, data = ng.fileio.bruker.read_pdata(path, all_components=True)            
        
        NS = self.acqus.NS
        RG = self.acqus.RG
        
        real = data[0]/(NS*RG)
        imag = data[1]/(NS*RG)
        self.espectro.set_real(real)
        self.espectro.set_imag(imag)
        self.espectro.set_size(real.shape)
        self.espectro.set_ppmAxis(self.crear_ppmAxis(self.procs))
        

    def crear_ppmAxis(self, procs):
        """
        """
        ppm = procs.ppm
        offset = procs.offset
        FTsize = procs.FTsize
        SpectralWidth = self.acqus.SW
        
        ultimo_ppm = offset - SpectralWidth
        ppmAxis = np.linspace(offset, ultimo_ppm, FTsize)
        return ppmAxis 
    
    
class DatosProcesados2D(DatosProcesados):

    def __init__(self, directorio, p_dir=1):
        DatosProcesados.__init__(self, directorio)
        self.proc2s = Procs(directorio, p_dir, dim2=True)
        self.espectro.set_ppmAxisInd(self.crear_ppmAxis(self.proc2s))
        self.espectro.crear_ppmGrid()
        
        
        