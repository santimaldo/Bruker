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



expnum = [3*i+1 for i in range(9)]
muestra = [str(exp) for exp in expnum]

rho = np.array([0,6,6,3,3,2,4,6,3])
phi = np.array([0,180,90,180,90,90,90,90,180])
phi_r = phi*np.pi/180
#------------------------------------------------------------
expnum = np.array([1,4,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91])

rho = np.array([0,6,  3, 2, 4, 6, 3,  4,  6, 3, 4, 2,6  ,3,  2,  4,  6,  3,  0,6  ,2  ,2 ,4  ,4 ,6 ,6  ,3 ,1,1 ])
phi = np.array([0,180,90,90,90,90,180,180,0, 0, 0, 0,270,270,180,180,180,180,0,180,270,90,270,90,90,180,90,0,90])
phi_r = phi*np.pi/180


expnum = np.array([31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,97,100,103,106,109,112,121,124])

rho = np.array([6, 3, 4, 2,6  ,3,  2,  4,  6,  3,  0,6  ,2  ,2 ,4  ,4 ,6 ,6  ,3 ,1,1 ,1  ,1  ,5,5 ,5  ,5  ,5  ,0,0])
phi = np.array([0, 0, 0, 0,270,270,180,180,180,180,0,180,270,90,270,90,90,180,90,0,90,180,270,0,90,180,270,270,0,0])
phi_r = phi*np.pi/180

# X:
r =      [-6 ,-6 ,-6 ,-6  ,-5 ,-4 ,-3 ,-2 ,-1 ,0 ,0  ,0  ,0  ,1 ,1  ,1  ,2 ,3 ,3  ,3  ,4 ,5  ,6 ]
expnum = [ 55, 64, 82, 133,106, 52, 58, 49, 94,61,121,124,154,88,142,145,40,34,148,151,37,100,31]

# Y:
# r =      [-6 ,-6 ,-6 ,-5  ,-4 ,-3 ,-2 ,-1 ,0 ,0  ,0  ,1 ,2 ,3 ,4 ,5  ,6 ]
# expnum = [ 55, 64, 82, 106, 52, 58, 49, 94,61,121,124,88,40,34,37,100,31]


P90_100_amp = []
P90_080_amp = []
P90_040_amp = []
P90_100_ph = []
P90_080_ph = []
P90_040_ph = []
# inicio un contador
cnt = 0
for expn in expnum:
  for jj in range(3):    
    directorio = path+str(expn+jj)+"/"
    datos = DatosProcesados(directorio)    
    NS = datos.acqus.NS  
    ph = datos.procs.phase
    
    spec = datos.espectro.real + 1j * datos.espectro.imag
    # spec = datos.espectro.real
    ppmAxis = datos.espectro.ppmAxis
    mod =np.abs(spec)
    baseline = np.mean(mod[:1000])
    mod = mod-baseline
    
    ini = 26
    fin = -16
    spec = spec[ppmAxis<ini]
    ppmAxis = ppmAxis[ppmAxis<ini]
    spec = spec[ppmAxis>fin]
    ppmAxis = ppmAxis[ppmAxis>fin]
    
    # plt.figure(expn)
    # plt.plot(ppmAxis, spec)
    # plt.xlim([20,-11])
    
    integral = np.trapz(mod)
    # integral = np.trapz(spec)
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
      print('graficando...')
      plt.figure(7896666)
      # plt.plot(ppmAxis, spec)
      # plt.plot(ppmAxis, ppmAxis*0,'k--')
      plt.plot(mod)
      plt.plot(mod*0,'k--')

      # plt.xlim([ini,fin])
    
    cnt += 1

#%%
r = np.array(r)
ind = np.where(r==0)

fases = [P90_100_ph, P90_080_ph, P90_040_ph]
fases_corr = []
for fase in fases:
  fase = np.array(fase)
  ph = np.mean(fase[ind])
  fases_corr.append(fase-ph)
  
integral = P90_100_amp 
integral = [i/integral[6] for i in integral]
fase = P90_100_ph
plt.figure(0)
plt.plot(r,integral,'o-')
plt.figure(10)
ph = [f-fase[6] for f in fase]
plt.plot(r,ph,'o-')


#%%
data = np.array([rho, phi, P90_100_amp, P90_080_amp, P90_040_amp, fases_corr[0], fases_corr[1], fases_corr[2]]).T
print(data)
# np.savetxt(savepath+'amplitudes.csv', data, delimiter=',', fmt='%i')
np.savetxt(savepath+'amplitudes.dat', data, fmt='%i')




fig = plt.figure(11123333)
area = 200
v = 90
#==============================================================================
ax = fig.add_subplot(231, projection='polar')
datos = P90_100_amp
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='hot_r', alpha=1, vmin=datos[0]*0.97) # , s=area
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ax = fig.add_subplot(234, projection='polar')
datos = np.array(P90_100_ph)-P90_100_ph[0]
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='seismic', alpha=1, vmin=-v, vmax=v) # , s=area
#==============================================================================
ax = fig.add_subplot(232, projection='polar')
datos = P90_080_amp
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='hot_r', alpha=1, vmin=datos[0]*0.97) # , s=area
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ax = fig.add_subplot(235, projection='polar')
datos = np.array(P90_080_ph)-P90_080_ph[0]
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='seismic', alpha=1, vmin=-v, vmax=v) # , s=area
#==============================================================================
ax = fig.add_subplot(233, projection='polar')
datos = P90_040_amp
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='hot_r', alpha=1, vmin=datos[0]*0.97) # , s=area
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ax = fig.add_subplot(236, projection='polar')
datos = np.array(P90_040_ph)-P90_040_ph[0]
c = ax.scatter(phi_r, rho, c=datos, s=area, cmap='seismic', alpha=1, vmin=-v, vmax=v) # , s=area

plt.show()


 #%%

P90_100_amp = np.array(P90_100_amp)


amp = P90_100_amp[phi%180==0] # la condicion es equivalente a phi==0 or phi == 180

Rho = rho[phi%180==0]
Phi = phi[phi%180==0]

x = Rho*np.exp(1j*Phi*np.pi/180.0)

plt.figure(909090909)
plt.plot(x,amp,'o')


