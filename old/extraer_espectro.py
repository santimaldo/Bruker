import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *

path  = "S:/Doctorado/Dendritas/116MHz/2019-09-11_Dendritas_SiO2/7/"
savepath = "S:/Doctorado/Dendritas/Analisis/2019-09_SiO2/"


path = '/media/santi/ESCARBATO/2020-10-05_Backup_lily/LiMetal/116MHz/2019-07-01_BB-Probe1.1-Test/100/'
savepath = "/home/santi/temp/"



path = 'S:/Doctorado/Carbones/300MHz/2021-08-30_CMK3_Activacion/{}/'
savepath = "S:/Doctorado/Carbones/analisis/2021-08_CMK3_Activacion/files/"
# filename = '2021-08-30_CMK3_Act_muestraD'


path = 'S:/Doctorado/Carbones/300MHz/2021-12-20_LiTFSI-TEGDME_CMK3/{}/'
savepath = "S:/Doctorado/Carbones/analisis/2021-12_CMK3/files/"



expnum = range(1,15)
expnum = range(2,10)
# expnum = [8]

for expn in expnum:
  print(expn)
  filename = str(expn)
  datos = DatosProcesados(path.format(expn))
  
  
  re = datos.espectro.real
  im = datos.espectro.imag
  
  
  ppmAxis = datos.espectro.ppmAxis
  re = re[np.abs(ppmAxis)<150]
  ppmAxis = ppmAxis[np.abs(ppmAxis)<150]
  
  
  dataexport = np.array([ppmAxis, re]).T
  np.savetxt(savepath+filename+'.dat', dataexport)
  
  print('graficando...')
  plt.figure(58)
  plt.plot(ppmAxis, re)
  #plt.plot(ppmAxis, im,'r')
  #plt.plot(mod, 'g')
plt.show()
