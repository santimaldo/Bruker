# coding=UTF-8

# librerias para la clase Datos
import numpy as np
from Acqus import *
from Procs import *
from PulseProg import *
# librerias para la clase DatosProcesados
from Espectro import *
from Fid import *
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from scipy.stats import linregress


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
        Es un objeto de la clase Acqus. Contiene los parametros con los cuales
        se realizo la medicion.
    procs : Procs
        Es un objeto de la clase Procs. Contiene los parametros con los cuales
        se procesaron los datos en TopSpin.
    pulseprog : PulseProg
        Es un objeto de la clase Procs. Contine informacion sobre el programa
        de pulso utilizado
    """

    def __init__(self, directorio, set_fid=False):

        self.directorio = directorio
        self.acqus = Acqus(directorio)
        self.pulseprog = self.getPulseProg()
        self.title = self.get_title()
        self.fid = Fid()

        self.nucleo = self.acqus.dic["NUC1"]
        self.gamma = self.set_gamma()  # rad/(s*T)

        # datos crudos
        if set_fid:
            self.set_fid()

    def UnMetodo(self):
        """
        @return  :
        @author
        """
        pass

    def get_title(self):
        """
        @return  :
        @author
        """
        title_file = "{}/pdata/1/title".format(self.directorio)
        with open(title_file, 'r') as file:
            data = file.read()
        return data

    def getPulseProg(self):
        """
        @return  : nombre del programa de pulsos
        """
        ppfile = f"{self.directorio}/pulseprogram"
        with open(ppfile, 'r') as file:
            # me quedo con la segunda linea del archivo
            data = file.readlines()[1]
            # eliminos los caracteres inicial y final: ; y \n
            data = data[1:-1]
        return data

    def set_gamma(self):
        """
        Establece el gamma de acuerdo al nucleo.
        unidades: rad/(s*T)

        fuente: https://en.wikipedia.org/wiki/Gyromagnetic_ratio
        """
        # paso a minuscula para evitar errores
        nucleo = self.nucleo.lower()
        if "h" in nucleo:
            gamma = 267.52218744
        elif "7li" in nucleo:
            gamma = 103.962
        elif "13c" in nucleo:
            gamma = 67.2828
        elif "19f" in nucleo:
            gamma = 251.815
        elif "23na" in nucleo:
            gamma = 70.761
        elif "31p" in nucleo:
            gamma = 108.291
        # -------------------
        gamma = gamma*1e6
        return gamma

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
            real = real[:, 69:]
            imag = imag[:, 69:]
            Npts = int(real[0, :].size)
        except IndexError:
            print('Es un espectro 1D')
            real = real[69:]
            imag = imag[69:]
            Npts = int(real.size)

        # defino el eje de tiempo
        DW = self.acqus.DW
        timeAxis = np.arange(Npts)*DW

        self.fid.set_timeAxis(timeAxis)
        self.fid.set_real(real)
        self.fid.set_imag(imag)
        self.fid.set_signal(real+1j*imag)
        self.fid.set_size(real.shape)
        self.fid.set_ppm(self.acqus.SFO1)

# ----------------------------HERENCIAS-----------------------------------------


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
        self.espectro.set_spec(real+1j*imag)
        self.espectro.set_size(real.shape)
        self.espectro.set_ppmAxis(self.crear_ppmAxis(self.procs, self.acqus))

    def crear_ppmAxis(self, procs, acqus):
        """
        """
        offset = procs.offset
        FTsize = procs.FTsize
        SpectralWidth = acqus.SW

        ultimo_ppm = offset - SpectralWidth
        ppmAxis = np.linspace(offset, ultimo_ppm, FTsize)
        return ppmAxis


class DatosProcesados2D(DatosProcesados):
    """
    espectroscopia 2D: T1, Difusion, EXSY, etc.

    signal: lista de integrales de los espectros (en la dmension directa)
    """

    def __init__(self, directorio, p_dir=1, ppmRange=None):
        DatosProcesados.__init__(self, directorio)
        self.proc2s = Procs(directorio, p_dir, dim2=True)
        self.acqu2s = Acqus(directorio, dim2=True)
        self.espectro.set_ppmAxisInd(
            self.crear_ppmAxis(self.proc2s, self.acqu2s))
        self.espectro.crear_ppmGrid()
        self.signal = None

        self.ppmRange = ppmRange

    def Integrar(self):
        """
        crea la lista de senales. Los datos ya tienen que estar procesados (xf2)

        Para ello utiliza una rango de ppms: 'ppmRange'. Si ninguno es provisto,
        utiliza todo el rango.
        """
        ppmRange = self.ppmRange

        ppmAxis = self.espectro.ppmAxis
        spec = self.espectro.real
        # defino rango de integracion
        if ppmRange is None:
            ppmRange = [self.espectro.ppmAxis[0], self.espectro.ppmAxis[-1]]

        midRange = (ppmRange[1]+ppmRange[0])/2
        semiRange = np.abs(ppmRange[1]-ppmRange[0])/2

        condicion = np.abs(ppmAxis-midRange) < semiRange
        newppm = ppmAxis[condicion]

        signal = []
        for n in range(spec[:, 0].size):
            spec1d = spec[n, :]
            spec1d = spec1d[condicion]
            # el signo menos de la integral es porque ppmAxis va
            # de mayor a menor.
            signal.append(-integrate.simps(spec1d, x=newppm))

        self.signal = np.array(signal)

        return 0


class DatosProcesadosT1(DatosProcesados2D):
    """
    NO ESTA TERMINADO

    atrubutos:
      tau:  lista de tiempos de T1 sat, en milisegundos
      signal:  S vs t procesado en topspin
    """

    def __init__(self, directorio, p_dir=1):
        DatosProcesados2D.__init__(self, directorio)

        self.tau = None
        self.signal = None
        self.factor_vd = None

        self.Crear_tau(directorio, p_dir)
        self.Integrar()

        self.tau_fit = None
        self.signal_fit = None
        self.T1params = None

    def Crear_tau(self, directorio, p_dir):
        """
        crea la lista de tau

        Busca el factor vd para multiplicar la lista vdlist
        """

        # Extraigo la lista de delays
        vdlist = np.loadtxt(f"{directorio}/vdlist")
        vdlist = vdlist*1000  # paso a ms

        # Extraigo el factor vd------------------------
        pulseprog = f"{directorio}/pulseprogram"
        data = []
        with open(pulseprog, 'rt') as f:
            for line in f:
                if line.lstrip().startswith('vd*'):
                    factor_vd = line.rstrip().split('*')
                    factor_vd = float(factor_vd[1])
                else:
                    continue
        self.factor_vd = factor_vd
        # ------------------------------------------------------------------------------

        self.tau = vdlist * factor_vd
        return 0

    def get_T1data(self, ppmRange=None):
        """
        Devuelve el par de arrays Tau y Signal.
        Debo darle el rango de integracion.
        Por defecto integra el rango completo.

        Atencion! Reescribe el atributo ppmRange
        """
        if ppmRange is not None:
            self.ppmRange = ppmRange
            self.Integrar()
        return self.tau, self.signal

    def T1fit(self):
        """
        Ajuste exponencial de T1

        solo implementado monoexponencial
        """
        def ExpDec1(tau, A, T1, y0, n=1):
            S = y0 + A * np.exp(-tau/T1)
            return S

        signal = self.signal
        tau = self.tau
        func = ExpDec1
        bounds = ([-np.infty, 0, -np.infty], [0, np.infty, np.infty])
        guess = (-signal[-1], self.factor_vd*1e3, signal[-1])  # A, T1, y0
        popt, pcov = curve_fit(func, tau, signal, guess, bounds=bounds)

        # calculo R^2
        residuals = signal - func(tau, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((-np.mean(signal))**2)
        r_squared = 1 - (ss_res / ss_tot)

        msg = f"Ajuste exponencial de T1:\n  " \
              f"T1 =  {popt[1]:.0f} ms\n  "    \
              f"Rsquared = {r_squared:.6f}"

        print(msg)

        self.tau_fit = np.linspace(0, self.tau[-1], 512)
        self.signal_fit = func(self.tau_fit, *popt)
        self.T1params = popt
        try:
            return self.tau_fit, self.signal_fit, residuals
        except ValueError:
            return self.tau_fit, self.signal_fit


# ------------------------------------------------------------------------------
class DatosProcesadosDiff(DatosProcesados2D):
    """
    NO ESTA TERMINADO. Por ahora funciona para secuencia STEbp

    atrubutos:
      bvalue:  lista de b value
      signal:  S vs b
      factor_b : factor de correccion del b value (calibracion de gradiente)
                 en unidades de 10^9 s^2/m

    """

    def __init__(self, directorio, p_dir=1, factor_b=1):
        DatosProcesados2D.__init__(self, directorio)

        self.bvalue = None
        self.signal = None
        self.factor_b = factor_b

        # parametros de la secuencia de gradiente
        self.delta = 2*self.acqus.P[30]*1e-3  # ms - para STE bipolar es 2*P30
        self.bigDelta = self.acqus.D[20]*1e3  # ms - para STE bipolar es D20
        self.Crear_bvalue()

        # quito el primer punto:
        Npts = self.acqu2s.TD
        self.espectro.IndirectDimensionSelect([1, Npts])
        self.bvalue = self.bvalue[1:]

        # Calculo los puntos S
        self.Integrar()

        self.signal_fit = None

    def Crear_bvalue(self):
        """
        crea la lista de bvalue
        """
        gpmax = self.acqus.gp
        Npts = self.acqu2s.TD
        delta = self.delta*1e-3  # s
        bigDelta = self.bigDelta  # s
        g0 = 12  # T/m
        gamma = self.gamma

        glist = np.linspace(0, g0*gpmax/100, Npts)

        # STE:
        bvalue = (gamma*glist*delta)**2 * (bigDelta-delta/3) * 1e-9
        self.bvalue = bvalue * self.factor_b
        return 0

    def get_Diffdata(self, ppmRange=None):
        """
        Devuelve el par de arrays Tau y Signal.
        Debo darle el rango de integracion.
        Por defecto integra el rango completo.

        Atencion! Reescribe el atributo ppmRange
        """
        if ppmRange is not None:
            self.ppmRange = ppmRange
            self.Integrar()
        return self.bvalue, self.signal

    def Diff1_fit(self):
        """
        Ajuste lineal de difusion
        """
        y = np.log(self.signal)
        x = self.bvalue

        slope, intercept, r, p, se = linregress(x, y)

        yfit = slope*x+intercept
        signal_fit = np.exp(yfit)
        residuals = self.signal - signal_fit

        D = -slope
        r_squared = r**2

        msg = f"Ajuste lineal de Difusion:\n \
          D =  {D:.3f} 10^-9 m^2/s\n \
          Rsquared = {r_squared:.6f}"

        print(msg)

        self.signal_fit = signal_fit
        try:
            return self.bvalue, self.signal_fit, residuals
        except ValueError:
            return self.bvalue, self.signal_fit

    # --------------------------------------------------------------Funciones
# %%
# plt.figure(1)
# plt.plot(tau, sig, 'o')
# plt.plot(t, s, '--')
# plt.show()
