import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Espectro import *
import scipy.integrate as integrate





path = 'S:/Doctorado/Carbones/analisis/2021-12_SAC/7T/2Dto1D/100ms/1/'
savepath = path
savefile = 'tmp.dat'

datos = DatosCrudos(path)
spec  = Espectro(fid=datos.fid)

# datos 2D

re = datos.fid.real
im = datos.fid.imag
t = datos.fid.timeAxis


re = spec.real
im = spec.imag
t = spec.ppmAxis


#s = (re + 1j* im) * np.exp(-20 * 1j * np.pi/180)
#re = np.real(s)
#im = np.imag(s)



plt.figure(0)
#plt.plot(re[15],'b-')
#plt.plot(im[15],'r-')
plt.plot(t,re,'b-')
plt.plot(t,im,'r-')

plt.show()





# savedata = np.array([t, re, im]).T
# np.savetxt(savepath+savefile,savedata)

#re = re[:,0:10]
#im = im[:,0:10]
#
#tp = np.linspace(1,300,300)
#
#integral_r = []
#integral_i = []
##for i in range(128):
##    Int_r = integrate.trapz(re[])
##    Int_i = integrate.trapz(im)
##    integral_r.append(Int_r)
##    integral_i.append(Int_i)
#
##plt.figure()
##plt.plot(tp, integral_r,'b-')
##plt.plot(tp, integral_i,'r-')
##plt.show()
#
#
#
#Int_r = integrate.simps(re)
#Int_i = integrate.simps(im)
#plt.figure(1)
#plt.plot(tp, Int_r,'b-')
#plt.plot(tp, Int_i,'r-')
#plt.show()
