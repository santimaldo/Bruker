# coding=UTF-8
import nmrglue as ng
import numpy as np

class Fid(object):

    """
    :version: 0
    :author: Santiago Agustin Maldonado Ochoa
    """

    """
    Clase Fid. Es el valor de la senal de RMN. Puede ser 1D o 2D.


    Parameteros
    -----------
    directorio : str
        Es la carpeta donde se encuentran los datos de Bruker.
        Los atributos son los parametros mas usados para el analisis de datos.

    Atributos
    ----------
    atributo : tipo
        Explicacion de que es.

    """
    def __init__(self):

        self.real = 0
        self.imag = 0
        self.signal = 0
        self.size = 1
        self.ppm = 0
        self.timeAxis = 0

    def set_real(self, real):
        self.real = real

    def set_imag(self, imag):
        self.imag = imag

    def set_signal(self, signal):
        self.signal = signal
        self.update_re_im()

    def set_size(self, size):
        self.size = size

    def set_ppm(self, ppm):
        self.ppm = ppm

    def set_timeAxis(self, timeAxis):
        self.timeAxis = timeAxis


    def em(self, em):
      """
      El famoso exponential multiplication de Bruker
      """
      t = self.timeAxis
      self.signal = self.signal*np.exp(-t*em)
      self.update_re_im()

    def RecortarTiempo(self, tmax):
      """
      tmax debe estar en segundos
      """

      self.signal = self.signal[self.timeAxis<tmax]
      self.timeAxis = self.timeAxis[self.timeAxis<tmax]
      self.update_re_im()


    def update_re_im(self):
      self.real = np.real(self.signal)
      self.imag = np.imag(self.signal)


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
