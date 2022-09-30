import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
# from VoigtFit import *

# directorio de datos
path = "S:/Doctorado/SaddleCoil/116MHz/2021-01-28_caracterizacion_saddle/"
# directorio de guradado
savepath = "S:/temp/"


# X:
r =      [-6 ,-6 ,-6 ,-5  ,-4 ,-3 ,-2 ,-1 ,0 ,0  ,0  ,1 ,2 ,3 ,4 ,5  ,6 ]
expnum = [ 55, 64, 82, 106, 52, 58, 49, 94,61,121,124,88,40,34,37,100,31]

# Y:
# r =      [-6 ,-6 ,-6 ,-5  ,-4 ,-3 ,-2 ,-1 ,0 ,0  ,0  ,1 ,2 ,3 ,4 ,5  ,6 ]
# expnum = [ 55, 64, 82, 106, 52, 58, 49, 94,61,121,124,88,40,34,37,100,31]



eje = []
amplitudes =[[],[],[]]
amplitudes_std = [[],[],[]]
fases = [[],[],[]]
fases_std = [[],[],[]]

P90_100_amp = []
P90_080_amp = []
P90_040_amp = []
P90_100_ph = []
P90_080_ph = []
P90_040_ph = []

r_temp = 666
# repeticiones:
amp_repeticiones = [[],[],[]]
ph_repeticiones = [[],[],[]]


for amp_rep in amp_repeticiones:
  amp_rep.append([])
for ph_rep in ph_repeticiones:
  ph_rep.append([])



# inicio un contador
cnt = 0  
for expn in expnum:      
  for jj in range(3):    
    directorio = path+str(expn+jj)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    ph = datos.procs.phase
    
    re = datos.espectro.real
    im = datos.espectro.imag    
    spec = re + 1j * im
    
    # spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    mod =np.abs(spec)
    baseline = np.mean(mod[:1000])
    mod = mod-baseline
    
    # ini = 26
    # fin = -16
    # spec = spec[ppmAxis<ini]
    # ppmAxis = ppmAxis[ppmAxis<ini]
    # spec = spec[ppmAxis>fin]
    # ppmAxis = ppmAxis[ppmAxis>fin]
    
    # plt.figure(expn)
    # plt.plot(ppmAxis, spec)
    # plt.xlim([20,-11])
    
    # integral = integrate.simps(mod)
    integral = np.abs(integrate.simps(spec))
    if ph<0:
      ph = ph+360
    # integral = np.trapz(spec)
    
    amp_repeticiones[jj][-1].append(integral)
    ph_repeticiones[jj][-1].append(ph)
    if jj==0:
      P90_100_amp.append(integral)
      P90_100_ph.append(ph)     
    elif jj==1:
      P90_080_amp.append(integral)
      P90_080_ph.append(ph)
    elif jj==2:
      P90_040_amp.append(integral)  
      P90_040_ph.append(ph)  
    
    # archivo_out = muestra[cnt]
    # dataexport = np.array([ppmAxis, spec]).T
    # np.savetxt(savepath+archivo_out+'.dat', dataexport)
    if expn>0:
      # print('graficando...')
      plt.figure(7896666)
      plt.plot(ppmAxis, np.real(spec))
      plt.plot(ppmAxis, np.imag(spec))

      # plt.xlim([ini,fin])
  r_temp = r[cnt-1]  
  if r[cnt]==r_temp:
    mismoRadio = True
    print(f"repito el {r[cnt]} ya que es igual a {r_temp}")
  else:
    print(f"-------------------------------r={r[cnt]}")
    eje.append(r[cnt])    
    for ii in range(3):
      amp_rep = amp_repeticiones[ii][-1]
      ph_rep = ph_repeticiones[ii][-1]
      if amp_rep!=[]:
        amp_rep=np.array(amp_rep)
        amplitudes[ii].append(np.mean(amp_rep))
        amplitudes_std[ii].append(np.std(amp_rep))
      if ph_rep!=[]:
        ph_rep=np.array(ph_rep)     
        fases[ii].append(np.mean(ph_rep))
        fases_std[ii].append(np.std(ph_rep))    
    # reinicio
    mismoRadio = False
    for amp_rep in amp_repeticiones:
      amp_rep.append([])
    for ph_rep in ph_repeticiones:
      ph_rep.append([])
  cnt += 1  

#%%
data = [r, P90_100_amp, P90_080_amp, P90_040_amp, P90_100_ph, P90_080_ph, P90_040_ph]
data = [np.array(dat) for dat in data]

for jj in range(-1,-4,-1):
  ind = np.where(data[0]==0) # busco donde r es igual a cero
  phase = data[jj]
  ph_centro = np.mean(phase[ind])
  print('fases: ', phase[ind],'|   media:', ph_centro)
  phase -= ph_centro
  data[jj]=phase
data = np.array(data).T
print(data)
# np.savetxt(savepath+'amplitudes.csv', data, delimiter=',', fmt='%i')
np.savetxt(savepath+'amplitudes.dat', data, fmt='%i')
#%%
 
colores = ['k','r','b'] 
for jj in range(3):
  plt.figure(1)
  amp = amplitudes[jj]
  amp_std = amplitudes_std[jj]
  plt.plot(eje, amp, 'o-', color=colores[jj])
  plt.errorbar(eje, amp, xerr=0.7, yerr=amp_std,color=colores[jj])

  # plt.figure(10)
  # ph = fases[jj]
  # ph_std = fases_std[jj]
  # plt.plot(eje, ph, 'o-', color=colores[jj])
  # plt.errorbar(eje, ph, xerr=0.7, yerr=ph_std,color=colores[jj])

