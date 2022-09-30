import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Exsy import *

from scipy.signal import savgol_filter


path  = "S:/Doctorado/Carbones/300MHz/2019-10-24_Carbones_MAS_EXSY/"
savepath = "S:/temp/"

# carbon + agua
mT =   [1 , 5 , 10, 25, 50, 100]
nexp = [25, 29, 26, 31, 27, 30 ]
picos = [1.9, -3.0]

guardar = False
# carbon  
# mT =   [1, 25, 50, 75, 100]
# nexp = [9, 15, 10, 13, 14]       
#picos = [0, -4.0]


# path  ="S:/Doctorado/Carbones/300MHz/2021-12-27_Carbones_MAS/"
# savepath = "S:/temp/"
# mT =   [10, 20, 50, 100, 200, 500, 1000]
# nexp = [25, 27, 24, 23 , 21 , 26 , 22]
# picos = [3.7, 2.7]


# path  ="S:/Doctorado/Carbones/300MHz/2022-05-06_Carbones_CMK3ACT/"
# savepath = "S:/temp/"
# mT = [100,10,1000,200,20,500,50,150,5,350,750,1,900,600,275,75,675,420,120,35] #ms
# nexp = np.arange(31,31+len(mT))
# picos = [3.6, 0.7]


integrales = []
espectros = []

for n in range(len(nexp)):
    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    datos = DatosProcesados2D(path+str(nexp[n])+"/")
    
    datos.espectro.ppmSelect2D([-8, 8])
    ##### anulo todos los valores menores a 0:
    datos.espectro.set_mask(0)
    
    exsy = Exsy(datos.espectro, 2)
       
    
    for i in range(len(picos)):
        for j in range(len(picos)):
            indices=(i,j)
            guess = (picos[i],picos[j])
            rango = (0.5,0.5) # rango de ppm en los cuales busca el maximo
            exsy.establecer_region(indices,guess, rango)
                  
    exsy.graficar_regiones(mT[n], Ncontour=100)
    integrales.append(exsy.integrar_regiones())
    
    if guardar:
        archivo_out = "mixTime_"+str(mT[n])+"ms"
        np.savetxt(savepath+archivo_out+'.dat', datos.espectro.real)
        np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
        np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)

    espectros.append(datos.espectro.real)
    X = exsy.ppmGridDir
    Y = exsy.ppmGridInd

#%%
    
# resta = espectros[1] - espectros[0]
# np.savetxt(savepath+'100ms-1ms.dat', resta)

i00 = []    
i10 = []
i01 = []
i11 = []
for n in range(len(nexp)):
    i00.append(integrales[n][0,0])
    i01.append(integrales[n][0,1])
    i10.append(integrales[n][1,0])
    i11.append(integrales[n][1,1])

i00 = np.array(i00)
i10 = np.array(i10)
i01 = np.array(i01)   
i11 = np.array(i11)

plt.figure(0)
#plt.plot(mT,i01,'o')
#plt.plot(mT,i10,'o')
plt.plot(mT,(i10+i01)/(i00+i11)/2,'o')
