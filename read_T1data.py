# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate




# directorio de datos
expn = 102
path = f"S:/Carbones/300MHz/2022-09-01_CarbonesCNEA_M4_LiTFSI/{expn}/"
# directorio de guradado
savepath= "S:/Carbones/analisis/2022-09_M4_LiTFSI-G2_100mM/filesT1/"
muestra = "tmp" #'1H_M4-NaOH' # "LiTFSI-G2_100mM_bulk", "M4-NaOH", "M4-HF"


# directorio de datos  ### M4 HF muestra SIN MOLER
# expn = 8
# path = f"S:/Carbones/300MHz/2018-09-13_Carbones_Glyma/{expn}/"
# savepath = "S:/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/T1/"
# muestra = "1H_M4-HF_SinMoler"

#### directorio de datos ### M4 HF muestra MOLIDA
expn = 11
path = f"S:/Carbones/300MHz/2018-11-28_Carbones_Glyma_MAS/{expn}/"
savepath = "S:/Carbones/analisis/2022-03_Glyme_AnalisisDatos2018/T1/"
muestra = "1H_M4-HF_Molida"



save = False
#rango de integracion
ppmRange = [20,-30]

datos = DatosProcesadosT1(path)
ppmAxis = datos.espectro.ppmAxis


spec = datos.espectro.real

re = datos.espectro.real[-1]
im = datos.espectro.imag[-1]


plt.figure(7532)
plt.plot(ppmAxis, re)
# plt.plot(ppmAxis, im)
plt.xlim(np.max(ppmAxis), np.min(ppmAxis))
plt.axhline(0, color='k')

#%%

tau, signal = datos.get_T1data(ppmRange)
tau_fit, signal_fit, residuals = datos.T1fit()



fig, axs = plt.subplots(2, 2)
fig.suptitle(muestra)
#-------------
axs[0, 0].plot(tau, signal, 'ko')
axs[0, 0].plot(tau_fit, signal_fit, 'r-')
text = f"$T_1 =$ {datos.T1params[1]:.0f} ms \n A = {datos.T1params[0]:.2f} \n $y_0 =$ {datos.T1params[2]:.2f}"
axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text,
          multialignment="left")
axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
#-------------
axs[1, 0].plot(tau, residuals, 'ko')
axs[1, 0].axhline(0, color='k', linestyle='--')
axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'Residuos')
#-------------
axs[0, 1].plot(tau, signal, 'ko')
axs[0, 1].plot(tau_fit, signal_fit, 'r-')
axs[0, 1].set(xlabel=r'$\tau$ [ms]')
axs[0, 1].set_xscale('log')
#-------------
axs[1, 1].plot(tau, residuals, 'ko')
axs[1, 1].axhline(0, color='k', linestyle='--')
axs[1, 1].set(xlabel=r'$\tau$ [ms]')
axs[1, 1].set_xscale('log')

for ax in axs.flat:
    ax.label_outer()

# guardo data:
if save:
  filename = f'{savepath}/{muestra}_T1.png'
  fig.savefig(filename)   # save the figure to file

  T1data = np.array([tau, signal]).T
  np.savetxt(f"{savepath}/{muestra}_T1.dat", T1data)

  data = np.array([ppmAxis, re, im]).T
  np.savetxt(f"{savepath}/{muestra}_ultimoEspectro.dat", data)







