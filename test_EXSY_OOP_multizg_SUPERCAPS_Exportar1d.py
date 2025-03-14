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
        i = (integrate.simpson(integrate.simpson(Matriz)))
    else:
        i = (integrate.simpson(integrate.simpson(Matriz, x=x), x=y))
    return i




# #############################################################################
# #############################################################################
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-03-10_insitu-LiTFSIaq-supercap/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\Supercaps\Analysis\2025-03_LiTFSI1M-aq_7Li-EXSY\exsy/-1V/"
muestra = f"Nexp{expn}"

mT = [0.1, 1, 10, 100, 100, 5, 50, 200, 20, 2, 80, 150]  # 1H
nexp = np.arange(80,91)


# reordeno - -------------------
# zipped_list = zip(mT[::2], nexp[::2])
zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
# sorted_list = sorted(zipped_list, reverse=True)
mT, nexp = np.array(sorted_list).T
# fin reordeno - ---------------
picos = [0, -5.9]
# parametrps -------------------
# semiancho de integracion (ppm)
semiancho = 0.5

# cmap = matplotlib.cm.Wistia
# cmap = matplotlib.cm.autumn
# cmap = matplotlib.cm.autumn_r
cmap = matplotlib.cm.YlOrBr_r
# nexp = [24]

filename = f"1H_EXSY_CMK3-ACT_semiancho{semiancho}"
rango = (6, -10)

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
    print(rf"{path}{int(nexp[n])}/", "  mT=", str(mT[n]))
    directorio = rf"{path}{int(nexp[n])}/"
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
                    filename = f"mT{mT[n]}"
                    np.savetxt(f"{savepath}{filename}.dat", data)
            else:
                idx_x = find_nearest(ppm_x, picos[j])
                spec1d = spec[:, idx_x]
                ppmaxis = ppm_x

    plt.legend()
    integrales.append(Int_n)
    # if mT[n]>=120: stop
