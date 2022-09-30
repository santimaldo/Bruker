import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
# from VoigtFit import *

# directorio de datos
path = "S:/Doctorado/LiMetal/116MHz/2021-04-06_SMC_Cogwheel/"
# path = "S:/Doctorado/LiMetal/116MHz/2021-02-02_SMC/"
# directorio de guradado


#abril
savepath = "S:/Doctorado/LiMetal/Analisis/2021-04_SMC/"
# expnum =  np.arange(3,11)
# k = np.arange(0.5,2.3,0.25)  # *1.16

expnum =  np.arange(20,27)
k = np.arange(0.5,2.1,0.25)


# febrero 
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N16_Delta1.5ms/"
# expnum =  np.arange(11,19)
# k = np.arange(0.5,2.3,0.25)

# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N16_Delta0.8ms/"
# expnum =  np.arange(19,27)
# k = np.arange(0.5,2.3,0.25)


# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N32_Delta1.5ms/"
# expnum =  np.arange(28,34)
# k = [0.5,1,1.25,1.5,1.75,2]


# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N16_Delta1.5ms-Muestra_ciclada/"
# expnum =  np.arange(41,48)
# k = np.arange(0.5,2.1,0.25)

# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N32_Delta1.5ms-Muestra_ciclada/"
# expnum =  np.arange(48,55)
# k = np.arange(0.5,2.1,0.25)

# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N64_Delta1.5ms-Muestra_Ciclada/"
# expnum =  np.arange(55,59)
# k = [1,1.1,1.2,1.3]

# #### zg pristino:
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/"
# expnum =  [10]
# k = [0]

# #### zg ciclado:
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/"
# expnum =  [40]
# k = [0]


cnt = 0
amplitud = []
ph = []
for expn in expnum:    
    directorio = path+str(expn)+"/"
    print(directorio)
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    phase = datos.procs.phase
    ph.append(phase)
    
    spec = datos.espectro.spec
    ppmAxis = datos.espectro.ppmAxis
    spec = spec[ppmAxis<300] 
    ppmAxis = ppmAxis[ppmAxis<300]
    spec = spec[ppmAxis>200] 
    ppmAxis = ppmAxis[ppmAxis>200]    
    
    archivo_out = 'k_{:.2f}'.format(k[cnt])
    dataexport = np.array([ppmAxis, np.real(spec), np.real(spec)/np.max(np.real(spec))]).T
    np.savetxt(savepath+archivo_out+'.dat', dataexport)
    
    dataexport = np.array([ppmAxis, np.abs(spec),np.abs(spec)/np.max(np.abs(spec))]).T
    # np.savetxt(savepath+archivo_out+'_mod.dat', dataexport)
    
    # amplitud.append(np.trapz(np.abs(spec)))
    amplitud.append(np.max(np.real(spec)))
    
    print('graficando...')
    plt.figure(0)
    # plt.plot(ppmAxis, spec)
    plt.plot(ppmAxis, np.real(spec)/np.max(np.real(spec)))
    
    cnt += 1


ax = plt.gca()
ax.invert_xaxis()
plt.show()


#%%
amplitud = np.array(amplitud)
ph=np.array(ph)

plt.figure(88888)
plt.plot(k,amplitud/np.max(amplitud),'o-')
plt.xlabel('k')

# plt.figure(88889)
# plt.plot(k,ph*np.pi/180,'ro-')
# plt.xlabel('k')