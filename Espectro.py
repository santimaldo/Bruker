# coding=UTF-8
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


class Espectro(object):
    """
    :version: 0
    :author: Santiago Agustin Maldonado Ochoa
    """

    """
  Clase Espectro. Es la transformada de fourier de la fid. Estos objetos
  pueden provenir de un procesamiento en TopSpin o bien de procesar una fid
  externamente


  Parameteros
  -----------
  directorio : str
      Es la carpeta donde se encuentran los datos de Bruker.

  Atributos
  ----------
  Los atributos son los parametros mas usados para el analisis de datos.

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

  Medotos
  -------

  abs():
        return: array
        devuelve el modulo del espectro.

  complex():
        return: array
        devuelve el expectro como numero complejo.

  ppmSelect(ppmRange):
        Redefine todos los parametros pero con el eje ppm recortado.

        parametros: list.
              ppmRange = [ppmInicial, ppmFinal]
  """

    def __init__(self, fid=None):  # , rango=None

        self.real = None
        self.imag = None
        self.spec = None
        self.ppmAxis = None
        self.ppmAxisInd = None
        self.ppmGridDir = None
        self.ppmGridInd = None
        self.size = None
        self.mask = None
        # self.rango=rango

        if fid is not None:
            self.CrearEspectro(fid)

        # print("RANGO::::", rango)
        # if self.rango is not None:
        #   if np.ndim(spec)==1:
        #     self.ppmSelect(rango)
        #   elif np.ndim(spec)==2:
        #     self.ppmSelect2D(rango)

    def set_real(self, real):
        self.real = real

    def set_imag(self, imag):
        self.imag = imag

    def set_spec(self, spec):
        self.spec = spec

    def set_size(self, size):
        self.size = size

    def set_ppmAxis(self, ppmAxis):
        self.ppmAxis = ppmAxis

    def set_ppmAxisInd(self, ppmAxisInd):
        self.ppmAxisInd = ppmAxisInd

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

    def ppmSelect(self, rango):
        """
        todavia no esta chequeado, no se si funciona
        """
        newppm = self.ppmAxis
        newppm = newppm[newppm <= max(rango)]
        newppm = newppm[newppm >= min(rango)]
        ini = np.where(self.ppmAxis == newppm[0])[0][0]
        fin = np.where(self.ppmAxis == newppm[-1])[0][0] + 1
        self.ppmAxis = newppm
        self.real = self.real[ini:fin]
        self.imag = self.imag[ini:fin]

    def ppmSelect2D(self, rango):
        newppm = self.ppmAxis
        newppm = newppm[newppm <= max(rango)]
        newppm = newppm[newppm >= min(rango)]
        dir_ini = np.where(self.ppmAxis == newppm[0])[0][0]
        dir_fin = np.where(self.ppmAxis == newppm[-1])[0][0] + 1
        self.ppmAxis = newppm

        newppm = self.ppmAxisInd
        newppm = newppm[newppm <= max(rango)]
        newppm = newppm[newppm >= min(rango)]
        ind_ini = np.where(self.ppmAxisInd == newppm[0])[0][0]
        ind_fin = np.where(self.ppmAxisInd == newppm[-1])[0][0] + 1
        self.ppmAxisInd = newppm

        self.real = self.real[ind_ini:ind_fin, dir_ini:dir_fin]
        self.imag = self.imag[ind_ini:ind_fin, dir_ini:dir_fin]

        self.crear_ppmGrid()

    def IndirectDimensionSelect(self, rango):
        """
        Selecciono los espectros en la dimension indirecta (xf2)

        rango : (int, int)
        """
        ini = int(rango[0])
        fin = int(rango[1])
        self.real = self.real[ini:fin, :]
        self.imag = self.imag[ini:fin, :]
        self.spec = self.real + 1j*self.imag
        self.crear_ppmGrid()

    def crear_ppmGrid(self):

        Dir = self.ppmAxis
        Ind = self.ppmAxisInd

        Dir, Ind = np.meshgrid(Dir, Ind)

        self.ppmGridDir = Dir
        self.ppmGridInd = Ind

    def set_mask(self, mask):
        """
        todos valores menores que la mascara, se fuerzan a cero.
        """
        self.mask = mask
        spec = self.real
        spec[spec <= mask] = 0.0
        self.real = spec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def CrearEspectro(self, fid, figure=None, return_angle=True):
        """
        Crear espectro a partir de una fid
        """

        # FID-----------------------------------------------------------------------
        ppm = fid.ppm  # Hz
        T2est = 0.14*1e-3  # T2est=0.14ms estimado con ancho de espectro. 2020-11-13
        t = fid.timeAxis
        signal = fid.signal
        dw = t[1]-t[0]
        NP = t.size

        # FOURIER-------------------------------------------------------------------
        freq = np.fft.fftshift(np.fft.fftfreq(NP, d=dw))
        ppmAxis = freq/ppm
        spec = np.fft.fftshift(np.fft.fft(signal))
        # corrijo la fase:
        spec, angle = autophase(ppmAxis, spec)

        self.ppmAxis = ppmAxis
        self.spec = spec
        self.real = np.real(spec)
        self.imag = np.imag(spec)

        if figure:
            plt.figure(figure)
            plt.subplot(1, 2, 1)
            plt.title('FID')
            plt.plot(t*1e3, np.real(fid)/np.max(np.real(fid)), linewidth=2)
            plt.xlabel('time [ms]')
            plt.yticks([])
            plt.subplot(1, 2, 2)
            plt.title('Spectrum')
            plt.plot(ppmAxis, np.real(spec)/np.max(np.real(spec)), linewidth=2)
            plt.xlabel('ppmAxis')
            plt.xlim([np.max(ppmAxis), np.min(ppmAxis)])
            plt.yticks([])

        if return_angle == True:
            return ppmAxis, spec, angle
        else:
            return ppmAxis, spec


def autophase(complex_data, x=None, method='minIntImag', precision=1):
    """
    Correccion automatica de fase, ya sea de espectro o de fid.

    complex_data es un espectro o fid, de la forma real + 1j* imag
    donde real, e imag son arrays del tipo float.

    si el metodo incluye la palabra 'fid', solo se considera el primer punto

    Usos:
    spec, angle = autophase(ppmAxis, spec)
    fid, angle = autophase(ppmAxis, fid, method='minIntImag')

    Metodos:
      + MinIntImagSpec (default)
        - espectros
        - minimiza el area de la parte imaginaria
      +

    """
    x_old = x
    complex_data_old = complex_data

    if x is None:
        # defino el eje x:
        x = np.linspace(complex_data.size)
    else:
        # ordeno el eje x:
        sort = np.argsort(x)
        x = x[sort]
        complex_data = complex_data[sort]

    angle = np.arange(0, 360, precision)

    SPECS = []
    IntImagSpec = []
    IntRealSpec = []

    if 'fid' in method.lower():
        Nmax = int(complex_data.size/(2**8))
        Nmax = 1
        x = x[:Nmax]
        complex_data = complex_data[:Nmax]

    for i in range(angle.size):
        Sp_try = complex_data*np.exp(-1j*angle[i]*np.pi/180)
        SPECS.append(Sp_try)
        if 'fid' in method.lower():
            intReal = np.real(Sp_try)[0]
            intImag = np.imag(Sp_try)[0]
        else:
            intReal = integrate.simps(np.real(Sp_try), x=x)
            intImag = integrate.simps(np.imag(Sp_try), x=x)

        IntRealSpec.append(intReal)
        IntImagSpec.append(intImag)
    IntImagSpec = np.array(IntImagSpec)
    IntRealSpec = np.array(IntRealSpec)

    if 'minintimag' in method.lower():
        # indice del minimo:
        idx = np.argmin(IntImagSpec)
        # ind_max = np.argmax(np.abs(np.real(complex_data)))
        # if complex_data[ind_max]<0:
        #   complex_data=-complex_data
    if 'maxintreal' in method.lower():
        # indice del maximo:
        idx = np.where(IntRealSpec == np.max(IntRealSpec))
        print(idx)
        idx = idx[0][0]
        ind_max = np.argmax(np.abs(np.real(complex_data)))
        # if complex_data[ind_max]<0:
        #   complex_data=-complex_data

    angle_opt = angle[idx]
    complex_data = complex_data_old * np.exp(-1j*angle_opt*np.pi/180)
    return complex_data, angle_opt
