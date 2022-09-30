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
    i = (integrate.simps(integrate.simps(Matriz, x=x),x=y))
  return i



###############################################################################
### SAC
path  ="S:/Doctorado/Carbones/300MHz/2021-12-23_Carbones_MAS/"
savepath = "S:/temp/"
savepath = "S:/Doctorado/Carbones/analisis/2022-05_Carbones_Sofi/CMK3_act_B_EXSY_1H/"
mT =   [1 , 50, 100, 250, 500, 750, 1000, 2000]
nexp = [18, 12, 13 , 14 , 15 , 16 , 17, 20]
# reordeno - -------------------
zipped_list = zip(mT, nexp)
sorted_list = sorted(zipped_list)
mT, nexp = np.array(sorted_list).T
# fin reordeno - ---------------
picos = [69,318]
# parametrps -------------------
# semiancho de integracion (ppm)
semiancho = 20

filename = f"1H_EXSY_SAC"


modulo=False
###############################################################################






# determino numero de filas y columnas del grafico-----------------------------

raiz = np.sqrt(len(nexp))
if int(raiz)**2 == len(nexp):
  filas = int(raiz)
  columnas = int(raiz)
else:
  filas = int(raiz//1) # parte entera
  columnas = int(raiz//1)+1
  
# otro caso limitado  
if filas*columnas<len(nexp):
  columnas+=1
#------------------------------------------------------------------------------

# chequeo que el semiancho no provoque solapamiento:---------------------------
ps = np.array(picos)
dists = np.abs(ps - np.roll(ps, -1))
if any(dists<2*semiancho):
  msj1=f"El semiancho {semiancho} ppm hace que se solapen algunos picos: {picos} ppm"
  msj2=f"\t El maximo semiancho posible es {np.min(dists/2)}"
  msj = f"{msj1}\n{msj2}"
  raise Exception(msj)
#------------------------------------------------------------------------------

fig_specs, ax_specs = plt.subplots()

integrales = []
espectros = []
for n in range(len(nexp)):
 

    print(path+str(nexp[n])+"/", "  mT=", str(mT[n]))
    directorio = f"{path}{nexp[n]}/"
    datos = DatosProcesados2D(directorio)
    # datos.espectro.ppmSelect2D(rango)
    #print(datos.title)

    if modulo:
      re = datos.espectro.real
      im = datos.espectro.imag    
      spec = np.abs(re+1j*im)
    else:
      spec = datos.espectro.real
    spec = spec[330:800,330:800]
    x = np.arange(spec[0,:].size)
    y = np.arange(spec[0,:].size)
    
    
    vmax=np.max(spec)
    plt.figure(151515)
    plt.title("mixing time: {} ms".format(mT[n])) 
    plt.subplot(filas, columnas,n+1)
    plt.contour(spec,30, colors=['k'], lw=0.1)
    
    spec1d = np.sum(np.real(spec), axis=0)
    ax_specs.plot(spec1d/np.max(spec1d), label=f"mixing time: {mT[n]:.0f} ms")

    Int_n = np.zeros([len(picos),len(picos)])        
    for i in range(len(picos)):       
      for j in range(len(picos)):
        xc = picos[i]           
        yc = picos[j]  
        
        xi= xc - semiancho
        xf= xc + semiancho
        yi= yc - semiancho
        yf= yc + semiancho
         
        spec_tmp = spec[yi:yf+1, xi:xf+1]
        x_tmp = x[xi:xf+1]
        y_tmp = y[yi:yf+1]
        
        plt.contour(x_tmp, y_tmp, spec_tmp, 5, cmap='inferno', lw=0.1)
        plt.xticks([])
        plt.yticks([])
        # plt.text(xc+semiancho, yc-semiancho, f'P{i+1}{j+1}',bbox=dict(facecolor='w', alpha=0.8))
        Int_n[i,j] = Integrar(spec_tmp)

    integrales.append(Int_n)
    # if mT[n]>=120: stop
           
integrales = np.array(integrales)
    # plt.figure(101010)
    # plt.subplot(filas, columnas,n+1)
    # plt.xticks([])
    # plt.yticks([])
    # plt.contour(spec,30, cmap='inferno')
    # plt.title("mixing time: {} ms".format(mT[n])) 
    # if guardar:
    #     archivo_out = "mixTime_"+str(mT[n])+"ms"
    #     np.savetxt(savepath+archivo_out+'.dat', datos.espectro.real)
    #     np.savetxt(savepath+archivo_out+'_ppmDir.dat', datos.espectro.ppmAxis)
    #     np.savetxt(savepath+archivo_out+'_ppmInd.dat', datos.espectro.ppmAxisInd)

    # espectros.append(datos.espectro.real)
    # X = exsy.ppmGridDir
    # Y = exsy.ppmGridInd

#%%

for i in range(len(picos)):
  for j in range(len(picos)):      
    if i<j:
      pii = integrales[:,i,i]
      pij = integrales[:,i,j]
      pji = integrales[:,j,i]
      pjj = integrales[:,j,j]
      # pii=1
      # pjj=1
      
      # cross peaks/ diag. peaks: medio, up y down
      cpu = pij/(pii+pjj)
      cpd = pji/(pii+pjj)
      cp = (pij+pji)/(pii+pjj)/2
      
      # CROSS PEAKS
      ii=i+1; jj=j+1
      plt.figure(ii*10+jj)
      plt.title(f"Picos {ii} ({picos[i]} ppm) y {jj} ({picos[j]} ppm)")
      plt.plot(mT, cpu, 'ro-', label= f"P{ii}{jj}/(P{ii}+P{jj})")
      plt.plot(mT, cpd, 'bo-', label= f"P{jj}{ii}/(P{ii}+P{jj})")
      plt.plot(mT, cp, 'ko-', label= f"[(P{ii}{jj}+P{jj}{ii})/2] / (P{ii}+P{jj})")
      # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
      plt.xlabel("Mixing time [ms]")
      plt.ylabel("Cross Peaks/Diag. Peaks [a.u]")
      plt.legend()
      
      plt.figure(ii*1000+jj)
      plt.title(f"Picos {ii} ({picos[i]} ppm) y {jj} ({picos[j]} ppm)")      
      plt.plot(mT, pij, 'o-', label= f"UP")
      plt.plot(mT, pji, 'o-', label= f"DOWN")
      # plt.plot(mT, (pii+pjj), 'ro-', label= f"(P{ii}+P{jj})")
      plt.xlabel("Mixing time [ms]")
      plt.ylabel("Cross Peaks")
      plt.legend()
      
      
      # guardado
      header = f"# Picos P{ii}:{picos[i]} ppm y P{jj} {picos[j]} ppm\n"
      # Up y Down son los que aparecen por arriba y por debajo de la diagonal
      header+= "# CrossPeaks/Diag.Peaks     : Pomedio de ambos picos cruzados\n"
      header+= "# CrossPeaksUp/Diag.Peaks   : Picos cruzados por arriba de la diagonal (Pij con i<j)\n"
      header+= "# CrossPeaksDown/Diag.Peaks : Picos cruzados por abajo de la diagonal  (Pij con i>j)\n"
      header+= f"# MixingTime_[ms]\tCrossPeaks/Diag.Peaks\tCrossPeaksUp/Diag.Peaks\tCrossPeaksDown/Diag.Peaks\tP{i}{j}(up)\tP{j}{i}(down)\tP{ii}\tP{jj}"
      data = np.array([mT, cp, cpu, cpd, pij, pji, pii, pjj]).T
      np.savetxt(f"{savepath}/{filename}_Picos{ii}{jj}.dat", data, header=header)

#%%      
# DIAGONAL PEAKS
for i in range(len(picos)):
  pii = integrales[:,i,i]      
  plt.figure(100)        
  plt.plot(mT, pii/pii[0], 'o-', label= f"P{i+1}: ({picos[i]} ppm)")
  plt.xlabel("Mixing time [ms]")
  plt.ylabel("Norm. Diag. Peaks")
  plt.legend()
        
#%%

plt.show()
# np.savetxt(f"{savepath}DATA2D_mT120.dat", spec)
# np.savetxt(f"{savepath}DATA2D_mT120_ppmAxis.dat", ppm_x)
# np.savetxt(f"{savepath}DATA2D_mT120_ppmAxisInd.dat", ppm_y)

# re = datos.espectro.real
# im = datos.espectro.imag    
# mod = np.abs(re+1j*im)
# np.savetxt(f"{savepath}DATA2D_mT120_ABS.dat", mod)