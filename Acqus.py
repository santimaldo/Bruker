# coding=UTF-8
import nmrglue as ng
import pathlib


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
    DW : float
        Dwell Time

  """

    def __init__(self, directorio, dim2=False):
        
        if not dim2:
            acq_file = 'acqus'
            self.acq_file = f'{directorio}/{acq_file}'
            self.dic = ng.bruker.read_acqus_file(directorio)[acq_file]
            self.NS = int(self.dic["NS"])
            self.D = self.dic["D"]
            self.D1 = self.D[1]
            self.P = self.dic["P"]
            self.P1 = self.P[1]
            self.PL1 = self.dic["PL"][1]
            self.RG = int(self.dic["RG"])        
            self.GPZ = self.dic["GPZ"]
            self.gp = self.GPZ[6]  # gradiente de la secuencia STE
            self.fecha = None
            self.hora = None
            self.getFechayHora()
        else:
            acq_file = 'acqu2s'
            self.acq_file = f'{directorio}/{acq_file}'
            self.dic = ng.bruker.read_acqus_file(directorio)[acq_file]
        # Here are the attributes that are common to both cases
        self.TD = int(self.dic["TD"])
        self.SW = self.dic["SW"]
        self.SFO1 = self.dic["SFO1"]
        sw_Hz = self.SW*self.SFO1
        self.DW = 1/(2*sw_Hz)  # en segundos
        

    def UnMetodo(self):
        """
        @return  :
        @author
        """
        pass

    def getFechayHora(self):
        FechayHora = pathlib.Path(self.acq_file).stat().st_mtime
        self.hora = FechayHora
