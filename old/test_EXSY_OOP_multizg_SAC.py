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
    i = (integrate.simps(integrate.simps(Matriz, x=y),x=x))
  return i


path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY/"
savepath = "S:/temp/"

# carbon + agua
mT =   [1 , 5 , 10, 25, 50, 100]
nexp = [25, 29, 26, 31, 27, 30 ]
picos = [1.9, -3.0]

guardar = False
# carbon  
mT =   [1, 25, 50, 75, 100]
nexp = [9, 15, 10, 13, 14]       
#picos = [0, -4.0]

#SAC
path  ="S:/Doctorado/Carbones/300MHz/2021-12-23_Carbones_MAS/"
savepath = "S:/temp/"
mT =   [1 , 50, 100, 250, 500, 750, 1000, 2000]
nexp = [18, 12, 13 , 14 , 15 , 16 , 17, 20]



#    y   x
# 1: 29,35  +- 5,10
# 2: 44,96  
#12: 44,35
#21: 29,96
xcentros=[35,96,35,96,35,96,35,96]
ycentros=[29,44,44,29,29,44,44,29]

integrales = []
espectros = []
for n in range(len(nexp)):

    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    datos = DatosProcesados2D(path+str(nexp[n])+"/")
    
    
    re = datos.espectro.real
    im = datos.espectro.imag
    
    spec = np.abs(re+1j*im)
    
    vcut = spec[650,650]*0.05 # recorto en el 5% del pico de microporos
    spec = spec[330:800,330:800]
    # en vez de recortar disminuyo exponencialmente para aquellos valores
    # por debajo de vcut
    # spec[spec<vcut]=np.exp(spec[spec<vcut]/vcut-1)*vcut

    # plt.figure(mT[n])
    # plt.contour(spec,5)
    # Int_n=[]
    # for i in range(len(xcentros)):
    #   xc = xcentros[i]
    #   yc = ycentros[i]

      
    #   sliced = spec[yc-3:yc+3, xc-10:xc+10]
    #   y_slice = np.arange(yc-3,yc+3)
    #   x_slice = np.arange(xc-10,xc+10)
    #   # plt.contourf(x_slice, y_slice, sliced,50)
    #   Int_n.append(Integrar(sliced))
    # integrales.append(Int_n)
    
    
    plt.figure(101010)
    plt.subplot(3,3,n+1)
    plt.xticks([])
    plt.yticks([])
    plt.contour(spec, 100, cmap='inferno')
    plt.title("mixing time: {} ms".format(mT[n])) 

    # if guardar:
    #     archivo_out = "mixTime_"+str(mT[n])+"ms"
    #     np.savetxt(savepath+archivo_out+'.dat', datos.espectro.real)
    #     np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
    #     np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)

    # espectros.append(datos.espectro.real)
    # X = exsy.ppmGridDir
    # Y = exsy.ppmGridInd
integrales = np.array(integrales)    

#%%

# p1, p2, p12, p21 = np.array(integrales).T
    

# plt.figure(0)
# plt.plot(mT,p1,'o')
# plt.plot(mT,p2,'o')
# plt.plot(mT,p12,'o')
# plt.plot(mT,p21,'o')
# # plt.plot(mT,(p1+p2)/2,'o')


# plt.figure(123456)
# plt.plot(mT,p12/(p1+p2),'o')
# plt.plot(mT,p21/(p1+p2),'o')
# plt.plot(mT,(p12+p21)/(p1+p2),'o')
# # plt.plot(mT,(p1+p2)/2,'o)

# #%%
# plt.figure(1023456)
# exsy = (p12+p21)/(p1+p2)
# exsy=exsy/np.max(exsy)
# error = np.abs((p12-p21)/(p1+p2))*2
# plt.plot(mT,exsy,'o')
# plt.errorbar(mT,exsy,yerr=error)

