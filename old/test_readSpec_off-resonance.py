import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from Datos import *
plt.rcParams.update({'font.size': 14})
# from VoigtFit import *

# directorio de datos
path = "S:/Doctorado/LiMetal/116MHz/2021-04-05_testBB/"
# directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-02_SMC/N32_Delta1.5ms/"




expnum = np.arange(221,242) # p1=21 us
nombre = r'$\tau_p = 21 \mu s$'
tp = 21e-6

expnum = np.arange(100,121) # p1=42 us
nombre = r'$\tau_p = 42 \mu s$'
tp = 42e-6

expnum = np.arange(200,221) # p1=84 us
nombre = r'$\tau_p = 84 \mu s$'
tp = 84e-6

expnum = np.arange(300,326) # p1=84 us
nombre = r'$\tau_p = 270 \mu s$'
tp = 270e-6

savepath = "S:/temp/"



cnt = 0
offset = []
amplitud = []
ph = []

ppm0 = 116.64171589
for expn in expnum:    
    directorio = path+str(expn)+"/"
    print(directorio)
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    phase = datos.procs.phase
    ph.append(phase)
    
    
    spec = datos.espectro.spec
    ppmAxis = datos.espectro.ppmAxis
    spec = spec[ppmAxis<10] 
    ppmAxis = ppmAxis[ppmAxis<10]
    spec = spec[ppmAxis>-10] 
    ppmAxis = ppmAxis[ppmAxis>-10]    
    
    # archivo_out = 'k_{:.2f}'.format(k[cnt])
    # dataexport = np.array([ppmAxis, np.real(spec), np.real(spec)/np.max(np.real(spec))]).T
    # np.savetxt(savepath+archivo_out+'.dat', dataexport)
    
    dataexport = np.array([ppmAxis, np.abs(spec),np.abs(spec)/np.max(np.abs(spec))]).T
    # np.savetxt(savepath+archivo_out+'_mod.dat', dataexport)
    
    offset.append((datos.acqus.SFO1-ppm0)*1e6/ppm0)
    amplitud.append(np.abs(integrate.simps(np.real(spec), ppmAxis)))
    # amplitud.append(np.abs(np.trapz(np.real(spec), x=ppmAxis)))
    
    # print('graficando...')
    # plt.figure(0)
    # # plt.plot(ppmAxis, spec)
    # plt.plot(ppmAxis, np.real(spec)/np.max(np.real(spec)))
    
    cnt += 1


# ax = plt.gca()
# ax.invert_xaxis()
# plt.show()


#%%
amplitud = np.array(amplitud)
offset = np.array(offset)


plt.figure(88888)
plt.plot(offset,amplitud/np.max(amplitud),'o--', label=nombre)
plt.xlabel(r'offset [ppm]')
plt.legend()

wnut= (np.pi/2) / tp

plt.figure(88889)
plt.plot(offset*ppm0/wnut, amplitud/np.max(amplitud),'o--', label=nombre)
plt.xlabel(r'$\Omega/\omega_{nut}$ [ppm]')
plt.legend()



# plt.figure(88889)
# plt.plot(k,ph*np.pi/180,'ro-')
# plt.xlabel('k')