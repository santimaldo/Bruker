import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Exsy import *
from scipy import integrate

from scipy.signal import savgol_filter


def Integrar(Matriz, x=None, y=None):
    """
    Integracion 2d
    """
    if x is None:
        i = (integrate.simps(integrate.simps(Matriz)))
    else:
        i = (integrate.simps(integrate.simps(Matriz, x=y), x=x))
    return i


path = "S:/Doctorado/Carbones/300MHz/2021-12-27_Carbones_MAS/"
path = "S:/NMRdata/2021_Carbones_Sofi/2021-12-27_Carbones_MAS/"  # oficin
savepath = "S:/tmp/"
mT = [100, 200, 500, 1000]
nexp = [32, 31, 30, 29]


#    y   x
# 1: 16,39  +- 3,9
# 2: 23,109
# 12: 23,39
# 21: 16,109
xcentros = [39, 109, 39, 109]
ycentros = [16, 23, 23, 16]

integrales = []
espectros = []


for n in range(len(nexp)):

    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    datos = DatosProcesados2D(path+str(nexp[n])+"/")

    spec = datos.espectro.real
    spec = spec[490:530, 920:1100]
    plt.figure(mT[n])
    plt.contour(spec, 5)
    Int_n = []
    for i in range(len(xcentros)):
        xc = xcentros[i]
        yc = ycentros[i]

        sliced = spec[yc-3:yc+3, xc-9:xc+9]
        y_slice = np.arange(yc-3, yc+3)
        x_slice = np.arange(xc-9, xc+9)
        plt.contourf(x_slice, y_slice, sliced, 50)

        Int_n.append(Integrar(sliced))
    integrales.append(Int_n)

    plt.figure(101010)
    plt.subplot(2, 2, n+1)
    plt.xticks([])
    plt.yticks([])
    plt.contour(spec[12:28, 16:150], 8, cmap='inferno')
    plt.title("mixing time: {} ms".format(mT[n]))

    # if guardar:
    #     archivo_out = "mixTime_"+str(mT[n])+"ms"
    #     np.savetxt(savepath+archivo_out+'.dat', datos.espectro.real)
    #     np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
    #     np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)

    # espectros.append(datos.espectro.real)
    # X = exsy.ppmGridDir
    # Y = exsy.ppmGridInd

# %%

p1, p2, p12, p21 = np.array(integrales).T


plt.figure(0)
plt.plot(mT, p1, 'o')
plt.plot(mT, p2, 'o')
plt.plot(mT, p12, 'o')
plt.plot(mT, p21, 'o')
# plt.plot(mT,(p1+p2)/2,'o')


plt.figure(123)
plt.plot(mT, p12/(p1+p2), 'o')
plt.plot(mT, p21/(p1+p2), 'o')
plt.plot(mT, (p12+p21)/(p1+p2), 'o')
# plt.plot(mT,(p1+p2)/2,'o)

# %%
plt.figure(1023)
exsy = (p12+p21)/(p1+p2)
exsy = exsy/np.max(exsy)
error = np.abs((p12-p21)/(p1+p2))*2
plt.plot(mT, exsy, 'o')
plt.errorbar(mT, exsy, yerr=error)
