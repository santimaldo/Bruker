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
from pathlib import Path

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

    def __init__(self, directorio, set_fid=False,
                 read_pp=True):

        self.directorio = directorio
        self.acqus = Acqus(directorio)
        self.pulseprog = None
        if read_pp:
            self.pulseprog = self.getPulseProg() 
        try:
            self.title = self.get_title()
        except FileNotFoundError:
            self.title = "No title file"
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
        ##### busco el pulseprogram--------------------------------
        # Buscar todos los archivos que contengan 'pulseprogram' en el nombre
        files = list(Path(self.directorio).glob("*pulseprogram*"))

        if not files:
            raise FileNotFoundError(f"No se encontró ningún archivo que contenga 'pulseprogram' en {self.directorio}")

        # Buscar coincidencia exacta con 'pulseprogram'
        exact_match = [f for f in files if f.name == "pulseprogram"]

        if len(exact_match) == 1:
            ppfile = str(exact_match[0])
        elif len(exact_match) > 1:
            print("Warning: Se encontraron múltiples archivos llamados exactamente 'pulseprogram'. Se usará el primero.")
            ppfile = str(exact_match[0])
        else:
            # Si no hay coincidencia exacta pero hay varios archivos
            if len(files) > 1:
                print("Warning: Se encontraron varios archivos que contienen 'pulseprogram', pero ninguno coincide exactamente.")
            ppfile = str(files[0])
        #----------------------------------------------------------

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
        elif "6li" in nucleo:
            gamma = 39.994
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
        elif "29si" in nucleo:
            gamma = -53.190
        elif "79br" in nucleo:
            gamma = 6.7256
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

        first_point = self.acqus.dic["GRPDLY"]+1
        # elimino los primeros puntos del tiempo muerto
        try:
            real = real[:, first_point:]
            imag = imag[:, first_point:]
            Npts = int(real[0, :].size)
        except IndexError:
            # print('Es un espectro 1D')
            real = real[first_point:]
            imag = imag[first_point:]
            Npts = int(real.size)
        # Normalizo por NS y RG
        NS = self.acqus.NS
        RG = self.acqus.RG
        real = real/(NS*RG)
        imag = imag/(NS*RG)

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

    def __init__(self, directorio, p_dir=1,
                 ppmRange=None,
                 read_pp=False,
                 nomrmalize_NS=True):
        Datos.__init__(self, directorio, read_pp=read_pp)
        self.p_dir = p_dir
        self.procs = Procs(directorio, p_dir)
        self.espectro = Espectro()
        self.ppmRange = ppmRange
        self.normalize_NS = nomrmalize_NS
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
        real = data[0]/(RG)
        imag = data[1]/(RG)
        if self.normalize_NS:
            real = real/(NS)
            imag = imag/(NS)    

        self.espectro.set_real(real)
        self.espectro.set_imag(imag)
        self.espectro.set_spec(real+1j*imag)
        self.espectro.set_size(real.shape)
        self.espectro.set_ppmAxis(self.crear_ppmAxis(self.procs, self.acqus))

    def crear_ppmAxis(self, procs, acqus, SF=None):
        """
        """
        ########## Implentacion original ################
        offset = procs.offset
        FTsize = procs.FTsize
        SpectralWidth = acqus.SW
        ultimo_ppm = offset - SpectralWidth
        ppmAxis = np.linspace(offset, ultimo_ppm, FTsize)
        ########## Otra implementacion ################
        # SWH = self.acqus.SWH
        # freq = np.linspace(SWH/2, -SWH/2, self.procs.dic["SI"])
        # ppmAxis = freq/self.acqus.SFO1
        # if SF is None:
        #     SF = self.procs.SF
        # ppmcenter = - 1e6*(self.acqus.SFO1 - SF)/self.acqus.SFO1
        # ppmAxis -= ppmcenter
        return ppmAxis
    
    def Integrar(self, ppmRange=None,absolute=False):
        """
        Calcula la señal integrada de un espectro 1D en el rango de ppm especificado.

        Si no se proporciona 'ppmRange', se utiliza el rango completo o el definido en el atributo de la clase.
        """
        ppmAxis = self.espectro.ppmAxis
        if absolute:
            spec = np.abs(self.espectro.spec)
        else:
            spec = self.espectro.real
        # defino rango de integracion
        if ppmRange is None:
            if self.ppmRange is not None:
                ppmRange = self.ppmRange
            else:
                print("WARNING: No se ha definido rango de integracion")
                ppmRange = [self.espectro.ppmAxis[0], self.espectro.ppmAxis[-1]]

        midRange = (ppmRange[1]+ppmRange[0])/2
        semiRange = abs(ppmRange[1]-ppmRange[0])/2

        condicion = abs(ppmAxis-midRange) < semiRange
        newppm = ppmAxis[condicion]

        signal = []
        
        spec1d = spec
        spec1d = spec1d[condicion]
        # el signo menos de la integral es porque ppmAxis va
        # de mayor a menor.
        signal = -integrate.simpson(spec1d, x=newppm)

        # Here I convert the list signal to a np array:
        signal = np.array(signal)
        self.signal = signal

        return signal


    def Mean_ppm(self, ppmRange=None,absolute=False):
        """
        Calcula el ppm promedio
        integrando de un espectro 1D en el rango de ppm especificado.

        mean_ppm = \int ppm * S(ppm) dppm / \int S(ppm) dppm
        Si no se proporciona 'ppmRange', se utiliza el rango completo o el definido en el atributo de la clase.
        """
        ppmAxis = self.espectro.ppmAxis
        if absolute:
            spec = np.abs(self.espectro.spec)
        else:
            spec = self.espectro.real
        # defino rango de integracion
        if ppmRange is None:
            if self.ppmRange is not None:
                ppmRange = self.ppmRange
            else:
                print("WARNING: No se ha definido rango de integracion")
                ppmRange = [self.espectro.ppmAxis[0], self.espectro.ppmAxis[-1]]

        midRange = (ppmRange[1]+ppmRange[0])/2
        semiRange = abs(ppmRange[1]-ppmRange[0])/2

        condicion = abs(ppmAxis-midRange) < semiRange
        newppm = ppmAxis[condicion]

        integral = []
        
        spec1d = spec
        spec1d = spec1d[condicion]
        # el signo menos de la integral es porque ppmAxis va
        # de mayor a menor.
        integral = integrate.simpson(newppm*spec1d, x=newppm)
        integral_total = integrate.simpson(spec1d, x=newppm) # pues no es una distribucion normalizada

        mean = integral / integral_total

        # Here I convert the list signal to a np array:
        mean = np.array(mean)


        return mean



class DatosProcesados2D(DatosProcesados):
    """
    espectroscopia 2D: T1, Difusion, EXSY, etc.

    signal: lista de integrales de los espectros (en la dmension directa)
    """

    def __init__(self, directorio, p_dir=1, ppmRange=None,
                 read_pp=False, extract_vdfactor=False):
        DatosProcesados.__init__(self, directorio, read_pp=read_pp)
        self.proc2s = Procs(directorio, p_dir, dim2=True)
        self.acqu2s = Acqus(directorio, dim2=True)
        self.espectro.set_ppmAxisInd(
            self.crear_ppmAxis(self.proc2s, self.acqu2s))
        self.espectro.crear_ppmGrid()
        self.signal = None
        self.directorio = directorio
        self.extract_vdfactor = extract_vdfactor

        self.ppmRange = ppmRange

    def Integrar(self, ppmRange=None,absolute=False):
        """
        crea la lista de senales. Los datos ya tienen que estar procesados (xf2)

        Para ello utiliza una rango de ppms: 'ppmRange'. Si ninguno es provisto,
        utiliza todo el rango.
        """
    

        ppmAxis = self.espectro.ppmAxis
        if absolute:
            spec = np.abs(self.espectro.spec)
        else:
            spec = self.espectro.real
        # defino rango de integracion
        if ppmRange is None:
            if self.ppmRange is not None:
                ppmRange = self.ppmRange
            else:
                print("WARNING: No se ha definido rango de integracion")
                ppmRange = [self.espectro.ppmAxis[0], self.espectro.ppmAxis[-1]]

        midRange = (ppmRange[1]+ppmRange[0])/2
        semiRange = abs(ppmRange[1]-ppmRange[0])/2

        condicion = abs(ppmAxis-midRange) < semiRange
        newppm = ppmAxis[condicion]

        signal = []
        for n in range(spec[:, 0].size):
            spec1d = spec[n, :]
            spec1d = spec1d[condicion]
            # el signo menos de la integral es porque ppmAxis va
            # de mayor a menor.
            signal.append(-integrate.simpson(spec1d, x=newppm))

        # Here I convert the list signal to a np array:
        signal = np.array(signal)
        self.signal = signal

        return signal
    
    def get_vdlist(self):
        """
        lee la lista de delays vdlist y los devuelve en ms
        """
        directorio = self.directorio
        # Extraigo la lista de delays
        if self.vdlist_path is None:
            vdlist = np.loadtxt(f"{directorio}/vdlist")
        else:
            vdlist = np.loadtxt(f"{directorio}/{self.vdlist_path}")
            print(f"vdlist leido de {directorio}/{self.vdlist_path}")
        vdlist = vdlist*1000 # paso a ms
        return vdlist


class DatosProcesadosT1(DatosProcesados2D):
    """
    NO ESTA TERMINADO

    atrubutos:
      tau:  lista de tiempos de T1 sat, en milisegundos
      signal:  S vs t procesado en topspin
    """

    def __init__(self, directorio, p_dir=1, vdlist_path=None):
        DatosProcesados2D.__init__(self, directorio)

        self.tau = None
        self.signal = None
        self.factor_vd = None
        self.vdlist_path = vdlist_path

        self.Crear_tau(directorio, p_dir)
        self.Integrar()

        self.tau_fit = None
        self.signal_fit = None
        self.T1params = None
        self.T1params_err = None

    def Crear_tau(self, directorio, p_dir):
        """
        crea la lista de tau

        Busca el factor vd para multiplicar la lista vdlist
        """

        # Extraigo la lista de delays
        vdlist = self.get_vdlist()

        # Extraigo el factor vd------------------------
        if self.extract_vdfactor:
            pulseprog = f"{directorio}/pulseprogram"
            data = []
            factor_vd = None
            with open(pulseprog, 'rt') as f:
                for line in f:
                    if line.lstrip().startswith('vd*'):
                        factor_vd = line.rstrip().split('*')
                        factor_vd = float(factor_vd[1])
                    else:
                        continue
            if factor_vd is None:
                print("WARNING: No se ha encontrado el factor vd. Seteando factor_vd = 1")
                factor_vd = 1
        else:
            factor_vd = 1
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
        size = min(self.tau.size, self.signal.size)
        self.tau = self.tau[:size]
        self.signal = self.signal[:size]
        return self.tau, self.signal

    def T1fit(self, model='mono'):
        """
        Ajuste exponencial de T1

        model : 'mono' o 'bi'
            - 'mono': y = y0 + A * exp(-tau/T1)
            - 'bi'  : y = y0 + A1 * exp(-tau/T11) + A2 * exp(-tau/T12)
        """
        signal = self.signal
        tau = self.tau

        if 'mono' in model.lower():
            def ExpDec(tau, A, T1, y0):
                return y0 + A * np.exp(-tau / T1)
            guess = (-signal[-1], self.factor_vd * 1e3, signal[-1])
            bounds = ([-np.inf, 0, -np.inf], [0, np.inf, np.inf])
            param_labels = ['A', 'T1', 'y0']

        elif 'bi' in model.lower():
            def ExpDec(tau, A1, T11, A2, T12, y0):
                return y0 + A1 * np.exp(-tau / T11) + A2 * np.exp(-tau / T12)
            guess = (-0.7*signal[-1], self.factor_vd * 0.5e3,
                    -0.3*signal[-1], self.factor_vd * 2e3,
                    signal[-1])
            bounds = (
                [-np.inf, 0, -np.inf, 0, -np.inf],
                [0, np.inf, 0, np.inf, np.inf]
            )
            param_labels = ['A1', 'T11', 'A2', 'T12', 'y0']
        else:
            raise ValueError("model debe ser 'mono' o 'bi'")

        # Ajuste
        popt, pcov = curve_fit(ExpDec, tau, signal, guess, bounds=bounds)

        # R^2
        residuals = signal - ExpDec(tau, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((signal - np.mean(signal))**2)
        r_squared = 1 - (ss_res / ss_tot)

        # === Salida y guardado ===
        self.tau_fit = np.logspace(np.log10(tau[0]),
                                np.log10(tau[-1]),
                                10*tau.size)
        self.signal_fit = ExpDec(self.tau_fit, *popt)
        self.T1params = popt
        self.T1stderr = np.sqrt(np.diag(pcov))
        self.T1params_err = 1.96 * self.T1stderr  # errores 95% confianza
        self.R2 = r_squared

        # Mensaje
        print(f"Ajuste {model}-exponencial de T1:")
        for name, val, err in zip(param_labels, popt, self.T1params_err):
            print(f"  {name} = {val:.2f} ± {err:.2f}")
        print(f"  R^2 = {r_squared:.6f}")

        return self.tau_fit, self.signal_fit, residuals


        

    def T2fit(self):
        """
        Ajuste exponencial de T2

        solo implementado monoexponencial
        """
        def ExpDec1(tau, A, T1, y0, n=1):
            S = y0 + A * np.exp(-tau/T1)
            return S

        signal = self.signal
        tau = self.tau
        func = ExpDec1
        bounds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])
        guess = (signal[0], self.factor_vd*1e3, 0)  # A, T1, y0
        popt, pcov = curve_fit(func, tau, signal, guess, bounds=bounds)

        # calculo R^2
        residuals = signal - func(tau, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((signal-np.mean(signal))**2)
        r_squared = 1 - (ss_res / ss_tot)

        msg = f"Ajuste exponencial de T2:\n  " \
              f"T2 =  {popt[1]:.0f} ms\n  "    \
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

    atributos:
      bvalue:  lista de b value
      signal:  S vs b
      factor_b : factor de correccion del b value (calibracion de gradiente)
                 en unidades de 10^9 s/m^2
      fit_results : [D, uD, r_squared] en unidades de 10^-9 m^2/s
      gpshape : str -- forma del pulso de gradientes. Opciones: 'rect', 'sin'
    """

    def __init__(self, directorio, p_dir=1, factor_b=1, bmax=np.inf, bmin=0, gpshape='sin'):
        # Herencia:
        DatosProcesados2D.__init__(self, directorio)

        # Inicializo atributos
        self.bvalue = None
        self.signal = None
        self.factor_b = factor_b
        self.bmax = bmax
        self.bmin = bmin
        self.absolute = None
        # parametros de la secuencia de gradiente
        # para STE bipolar, delta es 2*P30 y bigDelta es D20
        self.delta = 2*self.acqus.P[30]*1e-3  # ms - (originalmente en us)
        self.bigDelta = self.acqus.D[20]*1e3  # ms - (originalmente en  s)
        self.gpmax = self.acqus.GPZ[6]
        self.gpvalue = None

        self.Crear_bvalue(gpshape)
        # quito el primer punto:
        Npts = self.acqu2s.TD
        self.espectro.IndirectDimensionSelect([1, Npts])
        self.bvalue = self.bvalue[1:]

        # Calculo los puntos S
        self.Integrar()

        self.signal_fit = None

    def Crear_bvalue(self, gpshape):
        """
        crea la lista de bvalue
        """
        gpmax = self.acqus.gp
        Npts = self.acqu2s.TD
        delta = self.delta*1e-3  # s
        bigDelta = self.bigDelta*1e-3  # s
        g0 = 12  # T/m
        gamma = self.gamma
        
        # gp: gradiente porcentual; g: gradiente (T/m)
        gplist = np.linspace(0, gpmax, Npts)
        self.gpvalue = gplist[1:]
        glist = gplist/100*g0
        print(delta, bigDelta)

        # STE:
        if 'sin' in gpshape.lower():
            bvalue = (gamma*glist*delta/np.pi)**2 * (4*bigDelta-delta) * 1e-9
        elif 'rect' in gpshape.lower():
            bvalue = (gamma*glist*delta)**2 * (bigDelta-delta/3) * 1e-9

        self.bvalue = bvalue * self.factor_b
        return 0

    def get_Diffdata(self, ppmRange=None, absolute=False):
        """
        Devuelve el par de arrays Tau y Signal.
        Debo darle el rango de integracion.
        Por defecto integra el rango completo.

        Atencion! Reescribe el atributo ppmRange
        """
        self.absolute = absolute
        if ppmRange is not None:
            self.ppmRange = ppmRange
            self.Integrar(absolute=absolute)
        self.Recortar_datos()
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
        uD = se
        r_squared = r**2

        self.fit_results = [D, uD, r_squared]

        msg = f"Ajuste lineal de Difusion:\n \
          D =  {D:.8f} 10^-9 m^2/s\n \
          Rsquared = {r_squared:.6f}"

        print(msg)

        self.signal_fit = signal_fit
        try:
            return self.bvalue, self.signal_fit, residuals
        except ValueError:
            return self.bvalue, self.signal_fit

    def Recortar_datos(self):
        """
        Recorta los datos entre bmin y bmax.
        Por default, bmin=0 y bmax=inf.
        """
        if self.bmax == np.inf:
          bmax = self.bvalue[-1]
        else:
          bmax = self.bmax
        bmin = self.bmin
        if bmax != np.inf or bmin != 0:
            bvalue = self.bvalue
            signal = self.signal

            bmed = (bmax+bmin)/2.00
            semiancho = (bmax-bmin)/2.00
            condicion = np.abs(bvalue-bmed) <= semiancho
            # guardo en atributos:
            self.bvalue = bvalue[condicion]
            self.gpvalue = self.gpvalue[condicion]
            self.signal = signal[condicion]
        return 0
    # --------------------------------------------------------------Funciones
# %%
# plt.figure(1)
# plt.plot(tau, sig, 'o')
# plt.plot(t, s, '--')
# plt.show()
