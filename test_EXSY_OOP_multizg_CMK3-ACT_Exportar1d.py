import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Exsy import *
from Espectro import autophase
from scipy import integrate
import matplotlib.cm


from scipy.signal import savgol_filter


def Integrar(Matriz, x=None, y=None):
    """
    Integracion 2d
    """
    if x is None:
        i = (integrate.simps(integrate.simps(Matriz)))
    else:
        i = (integrate.simps(integrate.simps(Matriz, x=x), x=y))
    return i


# ###############################################################################
# path  ="S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY-reanalisis2022/"
# savepath = "S:/temp/"
# # mT =   [1,1,5,10,25,50,75,100] #ms
# # nexp = [8,9,5, 7,15,10,13, 14]
# mT =   [1,25,50,75,100] #ms
# nexp = [9,15,10,13, 14]
# # fin reordeno - ---------------
# picos = [-0.5,-4.3]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 1
# rango = (-12,6)
# ###############################################################################


##############################################################################
# CMK3ACTB Li
# path  ="S:/Doctorado/Carbones/300MHz/2022-05-19_Carbones_CMK3ACT/"
# path = "S:/NMRdata/2021_Carbones_Sofi/2022-05-19_Carbones_CMK3ACT/"  # compu Ofi
# savepath = path
# filename = "7Li_EXSY"
# mT = [1000, 300, 100, 600, 1, 25, 1250, 50,
#       450, 200, 800, 10, 150, 1500, 700, 375]
# nexp = np.arange(5, 5+2*len(mT), 2)
# # reordeno - -------------------
# zipped_list = zip(mT, nexp)
# sorted_list = sorted(zipped_list)
# mT, nexp = np.array(sorted_list).T
# # fin reordeno - ---------------
# picos = [-1.57, -3.75]
# # picos   = [-1.6,-2.8,-3.7]
# # parametrps -------------------
# # semiancho de integracion (ppm)
# semiancho = 0.5
# rango = (-7, 0)
# modulo = False
###############################################################################

# #############################################################################
# #############################################################################
# ###################CMK3ACTB 1H  # muestra "B". con cantidad de bulk correcta
# path = "S:/Doctorado/Carbones/300MHz/2022-05-12_Carbones_CMK3ACT/"  # Acer
path = "S:/Doctorado/Carbones/300MHz/2022-05-12_Carbones_CMK3ACT_reanalisis_09-2023/"  # Acer
# path = "S:/NMRdata/2021_Carbones_Sofi/2022-05-12_Carbones_CMK3ACT/"  # compu Ofi
# savepath = "S:/tmp"
savepath = "S:/tmp/"
# mT = [1, 100, 350, 1000, 10, 35, 600, 5, 20, 200, 75, 275,
#       50, 800, 150, 700, 900, 500, 420, 120, 1200, 1500]  # 1H
# nexp = np.arange(22, 22+2*len(mT), 2)

# para el slice a 0 ppm!!!:
mT = [1, 5, 10, 50, 100, 200, 500, 1000, 1500]  # 1H
nexp = [22, 36, 30, 46, 24, 40, 56, 28, 64]


# mT = [1, 5, 10, 50, 100]  # , 200, 500, 100, 1500]  # 1H
# nexp = [22, 36, 30, 46, 24]  # , 24, 40, 56, 28, 64]
# fases = [20, 20, 20, 5, 20]

# mT = [100, 200, 500]  # , 1000, 1500]  # 1H
# nexp = [24, 40, 56]  # , 28, 64]
# fases = [5, 10, -30]  # , 0, 0]


# reordeno - -------------------
# zipped_list = zip(mT[::2], nexp[::2])
zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
# sorted_list = sorted(zipped_list, reverse=True)
mT, nexp = np.array(sorted_list).T
# fin reordeno - ---------------
picos = [3.61, 0.63, -3.3]
# parametrps -------------------
# semiancho de integracion (ppm)
semiancho = 0.5

# cmap = matplotlib.cm.Wistia
# cmap = matplotlib.cm.autumn
# cmap = matplotlib.cm.autumn_r
cmap = matplotlib.cm.YlOrBr_r
# nexp = [24]

filename = f"1H_EXSY_CMK3-ACT_semiancho{semiancho}"
rango = (-20, 10)

modulo = False
##############################################################################


raiz = np.sqrt(len(nexp))
if int(raiz)**2 == len(nexp):
    filas = int(raiz)
    columnas = int(raiz)
else:
    filas = int(raiz//1)  # parte entera
    columnas = int(raiz//1)+1

# otro caso limitado
if filas*columnas < len(nexp):
    columnas += 1

# ------------------------------------------------------------------------------

# chequeo que el semiancho no provoque solapamiento:---------------------------
ps = np.array(picos)
dists = np.abs(ps - np.roll(ps, -1))
if any(dists < 2*semiancho):
    msj1 = f"El semiancho {semiancho} ppm hace que se solapen algunos picos: {picos} ppm"
    msj2 = f"\t El maximo semiancho posible es {np.min(dists/2)}"
    msj = f"{msj1}\n{msj2}"
    raise Exception(msj)
# ------------------------------------------------------------------------------


espectros1d = []
integrales = []
espectros = []
for n in range(len(nexp)):
    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    directorio = f"{path}{nexp[n]}/"
    datos = DatosProcesados2D(directorio)
    datos.espectro.ppmSelect2D(rango)
    # print(datos.title)

    if modulo:
        re = datos.espectro.real
        im = datos.espectro.imag
        spec = np.abs(re+1j*im)
    else:
        spec = datos.espectro.real
        im = datos.espectro.imag

    ppm_x = datos.espectro.ppmAxis
    ppm_y = datos.espectro.ppmAxisInd
    ppm_x_grid = datos.espectro.ppmGridDir
    ppm_y_grid = datos.espectro.ppmGridInd

    Int_n = np.zeros([len(picos), len(picos)])
    for i in range(len(picos)):
        for j in range(len(picos)):
            xc = picos[i]
            yc = picos[j]

            ppmx_tmp = ppm_x[np.abs(ppm_x-xc) < semiancho]
            ppmy_tmp = ppm_y[np.abs(ppm_y-yc) < semiancho]

            xi = np.where(ppm_x == ppmx_tmp[0])[0][0]
            xf = np.where(ppm_x == ppmx_tmp[-1])[0][0]
            yi = np.where(ppm_y == ppmy_tmp[0])[0][0]
            yf = np.where(ppm_y == ppmy_tmp[-1])[0][0]

            spec_tmp = spec[yi:yf+1, xi:xf+1]

            Int_n[i, j] = Integrar(spec_tmp, x=ppmx_tmp, y=ppmy_tmp)

            # Grafico espectros 1D
            color_index = n / (len(mT) - 1) * 0.9  # Normalized index
            color = cmap(color_index)

            if i == j:
                idx_y = find_nearest(ppm_y, picos[i])
                spec1d = spec[idx_y, :]
                ppmaxis = ppm_y
                if i == 0:
                    spec1dim = im[idx_y, :]
                    # s1d, _ = autophase(spec1d+1j*spec1dim, x=ppm_x)
                    # fase = fases[n]
                    # s1d = (spec1d+1j*spec1dim)*np.exp(-1j*fase*np.pi/180)
                    s1d = (spec1d+1j*spec1dim)
                    plt.figure(45138745953548)
                    plt.plot(ppm_x-3.6, s1d.real,
                             label=f"{mT[n]:.0f} ms")  # , color=color)
                    plt.xlim([5, -15])
                    plt.ylim([-200, 2000])
                    plt.xlabel(r"$\Delta\delta$ [ppm]")

                    data = np.array([ppm_x-3.6, spec1d]).T
                    filename = f"Chier_EXSY_mT{mT[n]}"
                    np.savetxt(f"{savepath}{filename}.dat", data)
            else:
                idx_x = find_nearest(ppm_x, picos[j])
                spec1d = spec[:, idx_x]
                ppmaxis = ppm_x

    plt.legend()
    integrales.append(Int_n)
    # if mT[n]>=120: stop
