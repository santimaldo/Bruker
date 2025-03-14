# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 18:58:34 2022

@author: santi
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import integrate

# Para graficar copiar estos parametros
# labels = {'xdata' : r'$x_{data}$ [units]',
#           'ydata' : r'$y_{data}$ [units]',
#           'xilt'  : r'$x_{ILT}$ [units]',
#           'yres'  : r'Residuals [units]'
#           'titulo':  'titulo'}


class ILT(object):
    """
    Transformada inveras de Laplace


    Parameteros
    -----------
    ydata    :: array
    xdata    :: array
    alpha    :: float
    kernel   :: function or str
    rango    :: 2-elements array
    Nilt     :: int
    xfit     :: array
    figure   :: int
    labels   :: dic
    savepath :: str   (debe terminar con / )
    """

    def __init__(self,
                 alpha=1e-3, kernel="T1sat", rango=(1, 1e4), Nilt=100,  # ILT
                 ydata=None, xdata=None,  # data
                 xfit=None,  # FIT
                 figure=None, labels=None,  # graficos
                 savepath=None, muestra=None):
        # Inicio listas de informacion:
        self.lista_de_muestras = []
        self.ii = 0  # aumenta cuando corro el DoTheStuff

        # directorio de guardado
        self.savepath = savepath

        # data
        self.ydata = ydata
        self.xdata = xdata
        self.muestra = muestra

        # ILT parameters
        self.alpha = alpha
        self.Nilt = Nilt
        self.rango = rango
        self.xilt = None
        self.yilt = None

        # Plot parameters
        self.figure = figure
        self.labels = labels

        # Kernel
        # Funciones disponibles para el kernel:
        self.kernel = kernel
        self.funciones = {"T1sat": self.T1sat}

        # Fit paremeters
        self.xfit = xfit
        # IMPLEMENTACION DE xfit:
        if xfit is not None:
            msg = "xfit no esta implementado todavia"
            raise Warning(msg)
        self.yfit = None
        self.residuals = None
        self.r_squared = None

        # graficos:
        if figure is not None:
            self.initialize_plot()


        if ydata is not None and xdata is not None:
            self.DoTheStuff(ydata, xdata, xfit=None, muestra=muestra)

    # ----------------------------------------------------------------------------

    def DoTheStuff(self, ydata, xdata, xfit=None, muestra=None, legend=False):
        print("#-------------------------------------")
        print("Calculando ILT ...")
        # data
        if muestra is None:
            self.muestra = f"Muestra_{self.ii}"
        else:
            self.muestra = muestra
        self.ydata = ydata
        self.xdata = xdata

        if xfit is not None:
            self.xfit = xfit
        else:
            if self.xfit is None:
                self.xfit = xdata

        # 1) Calcular
        self.Calculate()
        # 2) Graficar
        if self.figure is not None:
            self.plot()
            if legend:
                self.legend()
        # 3) Guardar
        if self.savepath is not None:
            self.save()

        # Actualizo informacion:
        self.ii += 1
        self.lista_de_muestras.append(muestra)

    def CrearKernel(self, kernel, Nk, Nilt):
        """
        Funcion para crear kernel, el cual es una matriz K de Nk filas y N columnas,
        donde:
            kernel:: puede ser dos tipos de inputs:
                      a) funcion : toma xdata y devuelve ydata
                      b) string  : secuencia tipica
                                    --> "T1sat", "T2", "STE", etc
            Nk    :: int  Numero de puntos de ydata, y del kernel
            Nilt  :: int  NUmero de puntos de la distribucion yilt

        Importante! la funcion kernel debe tomar dos variables, pera poder asociar
                    a xfata y xilt.
        """

        if isinstance(kernel, str):
            # Busca la funcion en el diccionario
            funcion = self.funciones[kernel]
            self.kernel = funcion
        else:
            funcion = kernel

        rango = self.rango
        xdata = self.xdata
        xilt = np.logspace(np.log10(min(rango)), np.log10(max(rango)), Nilt)

        K = np.zeros((Nk, Nilt))  # Npts_L: puntos de la transformada
        for ii in range(Nk):
            # K[ii,:] = 1 - np.exp(-tau[ii] / T)
            K[ii, :] = funcion(xdata[ii], xilt)

        self.xilt = xilt
        return K
    # ----------------------------------------------------------------------------

    def Calculate(self):
        '''
        Inversion de Laplace 1D
        '''
        ydata = self.ydata
        alpha = self.alpha
        Nilt = self.Nilt

        # Z : datos
        Z = ydata
        Z = np.reshape(Z, (len(Z), 1))
        # S : distribucion ILT
        S = np.zeros(Nilt)

        # creacion del Kernel
        Nk = Z.size
        K = self.CrearKernel(self.kernel, Nk, Nilt)

        KTK = K.T @ K  # @: matrix multiplication: @ = np.dot()
        KTZ = K.T @ Z
        ZZT = np.trace(Z @ Z.T)

        invL = 1 / (np.trace(KTK) + alpha)
        factor = 1 - alpha * invL

        Y = S
        tstep = 1
        lastRes = np.inf

        for iter in range(100000):
            term2 = KTZ - KTK @ Y
            Snew = factor * Y + invL * term2
            Snew[Snew < 0] = 0

            tnew = 0.5 * (1 + np.sqrt(1 + 4 * tstep**2))
            tRatio = (tstep - 1) / tnew
            Y = Snew + tRatio * (Snew - S)
            tstep = tnew
            S = Snew

            if iter % 1000 == 0:
                TikhTerm = alpha * np.linalg.norm(S)**2
                ObjFunc = ZZT - 2 * \
                    np.trace(S.T @ KTZ) + np.trace(S.T @ KTK @ S) + TikhTerm

                Res = np.abs(ObjFunc - lastRes) / ObjFunc
                lastRes = ObjFunc
                print(f'# It = {iter} >>> Residue = {Res:.6f}')

                if Res < 1E-5:
                    break

        self.yilt = S[:, 0]
        self.fit()

        return self.xilt, self.yilt
    
    def initialize_plot(self):
        """
        Metodo para inicializar el grafico
        """
        fig = plt.figure(num=self.figure, figsize=(12, 8))
        ax0 = plt.subplot2grid((3, 2), (0, 0), rowspan=2)
        ax1 = plt.subplot2grid((3, 2), (2, 0))
        ax2 = plt.subplot2grid((3, 2), (0, 1), rowspan=3)
        self.fig = fig
        self.axes = [ax0, ax1, ax2]
        return fig, [ax0, ax1, ax2]

    def fit(self):
        '''
        Ajuste del decaimiento a partir de la distribucion
        '''
        xdata = self.xdata
        ydata = self.ydata
        xilt = self.xilt
        yilt = self.yilt

        xfit = self.xfit
        M = []
        # IMPLEMENTACION DE xfit: aca la suma deberia ser range(xfit.size)
        for i in range(xfit.size):
            m = 0
            for j in range(xilt.size):
                m += yilt[j] * self.kernel(xfit[i], xilt[j])
            M.append(m)
        yfit = np.array(M)
        self.yfit = yfit

        # calculo R^2
        residuals = ydata - yfit  # IMPLEMENTACION DE xfit: ojo con las dims
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        r_squared = 1 - (ss_res / ss_tot)

        self.residuals = residuals
        self.r_squared = r_squared

        return yfit

    # -----------------------------Grafico--------------------------------------

    def plot(self):
        """
        Metodo para plotear la transformada
        """

        plt.rcParams.update({'font.size': 16})

        if self.labels is None:
            labels = {'xdata': r'$x_{data}$ [units]',
                      'ydata': r'$y_{data}$ [units]',
                      'xilt': r'$x_{ILT}$ [units]',
                      'yres': r'Residuals [units]',
                      'titulo': None}
        else:
            labels = self.labels

        xdata = self.xdata
        ydata = self.ydata
        xilt = self.xilt
        yilt = self.yilt
        xfit = self.xfit
        yfit = self.yfit
        residuals = self.residuals
        r_squared = self.r_squared

        lw = 3
        fig = self.fig
        ax0, ax1, ax2 = self.axes

        if labels['titulo'] is not None:
            fig.suptitle(labels['titulo'])
        else:
            fig.suptitle(self.muestra)

        ax0.plot(xdata, ydata, 'o')
        ax0.plot(xfit, yfit, 'r-', lw=lw)
        ax0.set_xlabel(labels['xdata'])
        ax0.set_ylabel(labels['ydata'])

        ax1.plot(xfit, residuals, 'o')
        ax1.axhline(0, color='k')
        ax1.set_xlabel(labels['xdata'])
        ax1.set_ylabel(labels['yres'])
        textstr = fr"$R^{{2}} = {r_squared:5f}$"
        ax1.text(0.6, 0.90, textstr, transform=ax1.transAxes, fontsize=14,
                 verticalalignment='top', bbox={'alpha': 1, 'facecolor': 'white'})

        cumulative = integrate.cumulative_trapezoid(yilt[2:-1])
        cumulative = cumulative / cumulative.max()

        ax2.semilogx(xilt[3:-1], cumulative, 'k', lw=lw/2)
        ax2.semilogx(xilt, yilt/yilt.max(), lw=lw, label=self.muestra)
        ax2.set_xlabel(labels['xilt'])
        ax2.yaxis.tick_right()
        ax2.set_yticks(np.linspace(0, 1, 11))
        ax2.grid(axis='y', ls='--')

        fig.tight_layout()
        axs = [ax0, ax1, ax2]

        return fig, axs

    def legend(self):
        """
        Metodo agregar legend si tengo mas de una curva
        """
        ax2 = self.axes[2]
        # ax2.legend(self.lista_de_muestras)
        ax2.legend()
    # ---------------------------------------guardado---------------------------

    def save(self):
        if self.figure is not None:
            self.fig.savefig(f"{self.savepath}{self.muestra}.png")

        data = np.array([self.xilt, self.yilt]).T
        header = f"ILT Distribution\n"\
                 f"Kernel\t:\t{self.kernel.__name__} \n"\
                 f"Xilt\t Yilt\t"
        filename = f'{self.savepath}{self.muestra}_ILT-Dist.dat'
        np.savetxt(filename, data, header=header)
        # - - - - - - - - - - - -
        data = np.array([self.xdata, self.ydata,
                         self.xfit, self.yfit,
                         self.residuals]).T
        header = f"Data and ILT fit\n"\
                 f"Kernel\t:\t{self.kernel.__name__} \n"\
                 f"R^2\t:\t{self.r_squared:6f}\n"\
                 f"Xdata\t Ydata\t Xfit\t Yfit\t residuals"
        filename = f'{self.savepath}{self.muestra}_ILT-fit.dat'
        np.savetxt(filename, data, header=header)

    def ILTFT(self, ydata2d, xdata, Nsegments=32):

        # First, we determine the length of the
        # segments in the direct dimention
        # we will use to calculate the ILT
        Nind, Ndir = ydata2d.shape  # number of data points
        segm_len = Ndir // Nsegments
        Nilt = self.Nilt
        iltft = np.zeros([Nilt, Nsegments])
        ydata2d_reduced = np.zeros([Nind, Nsegments])
        # set None to prevent multiple plots
        self.figure = None
        
        # Compute the number of columns per segment
        columns_per_segment = Ndir // Nsegments  # Integer division
        # Trim extra columns if needed
        trimmed_columns = columns_per_segment * Nsegments
        ydata2d_trimmed = ydata2d[:, :trimmed_columns]  # Remove excess columns
        # Reshape and compute the mean over the last dimension
        ydata2d_reduced = ydata2d_trimmed.reshape(Nind, Nsegments, columns_per_segment).mean(axis=2)


        for nseg in range(Nsegments):
            print(f"Processing segment {nseg+1} of {Nsegments}")
            ydata = ydata2d_reduced[:,nseg]
            self.DoTheStuff(ydata, xdata)
            iltft[:,nseg] = self.yilt

        return iltft, ydata2d_reduced

    ##################### Funciones Kernel ####################################
    # xdata es el eje temporal de la senal.
    # xilt es la variable de la ILT, es decir, T1, T2, D o lo que corresponda.

    def T1sat(self, xdata, xilt):
        """
        Saturation Recovery
        """
        return 1 - np.exp(-xdata / xilt)


# %%---------------------------------------------------------------------------
if __name__ == "__main__":

    # datos inventados:
    a1 = np.random.rand()
    a2 = np.random.rand()
    t11, t12, t21, t22 = np.random.rand(4)*1000 + 100
    xdata = np.linspace(0, 5000, 1024)
    ydata1 = a1 * np.exp(-xdata/t11) + (1-a1) * \
        np.exp(-xdata/t12) + 0.05*np.random.rand(1024)
    ydata2 = a2 * np.exp(-xdata/t21) + (1-a2) * \
        np.exp(-xdata/t22) + 0.05*np.random.rand(1024)

    Npts_L = 256
    alpha = 1e-3

    # funcion de kernel
    def T2(x, tt):
        return np.exp(-x/tt)

    # Inicializo la clase ILT
    ilt = ILT(alpha=alpha, rango=(1e1, 1e5), kernel=T2, Nilt=100,
              figure=2, savepath='c:/tmp/')
    # calculo la ILT para el conjunto de datos 1
    ilt.DoTheStuff(ydata1, xdata)
    # calculo la ILT para el conjunto de datos 2
    ilt.DoTheStuff(ydata2, xdata, muestra="data_2")
    ilt.legend()

    xilt = ilt.xilt
    yilt = ilt.yilt
    xfit = ilt.xfit
    yfit = ilt.yfit
    residuals = ilt.residuals
    r_squared = ilt.r_squared
    #%%
    data = np.loadtxt(r"c:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir\Nexp_72.dat")
    xdata = data[:, 0]
    ydata1 = -(data[:, 1]-1)

    data = np.loadtxt(r"c:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir\Nexp_73.dat")
    xdata = data[:, 0]
    ydata2 = -(data[:, 1]-1)

    data = np.loadtxt(r"c:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir\Nexp_74.dat")
    xdata = data[:, 0]
    ydata3 = -(data[:, 1]-1)

    data = np.loadtxt(r"c:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\ir\Nexp_75.dat")
    xdata = data[:, 0]
    ydata4 = -(data[:, 1]-1)
    # funcion de kernel
    def T2(x, tt):
        return np.exp(-x/tt)
    
    Npts_L = 256
    alpha = 1e-2

    # Inicializo la clase ILT
    ilt = ILT(alpha=alpha, rango=(1e-2, 1e1), kernel=T2, Nilt=100,
              figure=2, savepath='c:/tmp/')
    # calculo la ILT para el conjunto de datos 1
    ilt.DoTheStuff(ydata1, xdata, muestra="-1V")
    # calculo la ILT para el conjunto de datos 2
    ilt.DoTheStuff(ydata2, xdata, muestra="1V")
    ilt.DoTheStuff(ydata3, xdata, muestra="0V")
    ilt.DoTheStuff(ydata4, xdata, muestra="-1V again")
    ilt.legend()
# %%
