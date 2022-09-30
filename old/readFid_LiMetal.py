import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from Datos import *
from VoigtFit import *




# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2021-11-08_cabezal_litio_senal_vs_z/"


# celda S1
alturas = [0 , 0.75, 1.5, 3 , 5 ] 
expnum  = [21, 26  , 24 , 25, 22 ]
# directorio de guradado
savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_B1_vs_z/Swagelok1/"

# # celda S2
# alturas = [0, 0.75, 1.5, 3, 4.5 ] 
# expnum  = [5, 6   , 4  , 2, 3 ]
# # directorio de guradado
# savepath = "S:/Doctorado/LiMetal/Analisis/2021-11_B1_vs_z/Swagelok2/"

# Aluminas
# directorio de datos
path  = "S:/Doctorado/LiMetal/116MHz/2021-11-23_Aluminas/"
expnum  = [4]
filenames = ["LiMetal"]
# directorio de guradado
savepath = "S:/temp/"

# correccion de fase, grados
phase = 10

# t2 estrella (s)
T2est = 0.14e-3


# plt.figure(123568)
# plt.xlim([0,0.001])
# plt.hlines(0,0,0.001)


integrales = []
integrales_err = []
jj=0
for expn in expnum:
    directorio = path+str(expn)+"/"
    datos = DatosCrudos(directorio)
    
    timeAxis = datos.fid.timeAxis
    re       = datos.fid.real
    im       = datos.fid.imag
    
   
    s = re + 1j* im 
    
    s_new = s*np.exp(1j*phase*np.pi/180)
    
    re_new = np.real(s_new)
    im_new = np.imag(s_new)
   
    # print('graficando...')    
    # plt.plot(timeAxis, im_new, 'r')
    # plt.plot(timeAxis, re_new, 'k')
    
    jj+=1


S_abs = np.abs(s)
maxx = np.mean(S_abs[0:4])
S_abs=S_abs/maxx
base = np.mean(S_abs[-50:])
S_abs -= base



plt.figure(21)
plt.plot(timeAxis, S_abs, 'o')
plt.plot(timeAxis, np.exp(-timeAxis/T2est))