# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi
"""

import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Autophase import autophase
import scipy.integrate as integrate
import matplotlib.ticker as ticker
from Laplace import ILT
from scipy.optimize import curve_fit

save = False
expns = [9, 22]
samples = ["19F_LP57-neat", "19F_LP57-Mn(TFSI)2-8mM"]
substract_offset = False

expns = [14, 25]
samples = ["7Li_LP57-neat", "7Li_LP57-Mn(TFSI)2-8mM"]
substract_offset = True
#-------------------- directorios
path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-10-13_DRinsitu_LP57-MnTFSI/"
# directorio de guradado
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\TMions\Analysis\2025-10_in-situ_relaxation/"


cpmgs = []
times = []
for i, [expn, sample] in enumerate(zip(expns, samples)):
    # --------------------------- Extraigo datos
    datos = Datos(f"{path}{expn}", set_fid=True)
    re = datos.fid.real
    im = datos.fid.imag
    echotime = datos.acqus.dic["D"][20] * 2 # *2 es porque D20 es el tiempo entre el pulso y el eco
    # i0 = datos.acqus.dic["GRPDLY"] + 1  # primer punto valido del FID
    timeAxis = np.arange(1, re.size+1) * (echotime) 
    fid = re+1j*im
    mod = np.abs(fid)
    # angle= -61
    # fid = fid * np.exp(1j*angle*np.pi/180)
    fid, angle = autophase(timeAxis, fid, method='minIntImag')



    re = fid.real
    im = fid.imag
    # -----------------------------------------------
    cpmgs.append(re)
    times.append(timeAxis)
    #%%
    titulo = f"sample: {sample}"

    fig1d, ax1d = plt.subplots(1, 1)
    ax1d.axhline(0, color='gray')
    ax1d.plot(timeAxis*1000, re, 'ko-', lw=2, label='Real')
    ax1d.plot(timeAxis*1000, im, 'ro-', lw=2, label='Imag')
    ax1d.plot(timeAxis*1000, mod, 'grey', lw=2, alpha=0.4, label='Abs')
    # ax1d.set_xlim(0,0.05)
    # ax1d.plot(timeAxis, np.real(fid), 'k--', lw=2, label='Real')
    # ax1d.plot(timeAxis, np.imag(fid), 'r--', lw=2, label='Imag')

    ax1d.set_xlabel("Time [ms]")
    ax1d.set_title(titulo)
    ax1d.legend()
    plt.axhline(0, color='k')
# -----------------------------------------------
#%% ILT
Npts_L = 100
alpha = 1e-2

# funcion de kernel
def T2(x, tt):
    return np.exp(-x/tt)

labels = {'xdata': r'Time [s]',
          'ydata': r'Signal amplitude [a.u]',
          'xilt': r'$T_{2}$ [s]',
          'yres': r'Residuals [a.u]',
          'titulo': "CPMG"}

# Inicializo la clase ILT
ilt = ILT(alpha=alpha, rango=(1e-4, 1e1), 
          kernel=T2, Nilt=Npts_L,
          figure=2, savepath=savepath,
          labels=labels,
          yscale='log'
          )
# calculo la ILT para el conjunto de datos 1
for ii in range(len(cpmgs)):
    ydata = cpmgs[ii]/cpmgs[ii][0]
    if substract_offset:
        ydata = ydata - np.mean(ydata[-100:])  # offset
    ilt.DoTheStuff(ydata, times[ii], muestra=samples[ii])
ilt.legend()
#%% biexponential fit

# Obtener la paleta Set1
cmap = plt.get_cmap("Set1")
colors = [cmap(i) for i in range(cmap.N)]  # 9 colores

# Función biexponencial con offset
def biexp_offset(t, A1, T1, A2, T2, y0):
    return A1 * np.exp(-t/T1) + A2 * np.exp(-t/T2) + y0

# Listas para resultados
T1short_list = []
T1long_list = []
P1short_list = []
offset_list = []

font_size = 10  # tamaño de la fuente

# Crear figura con dos subplots (curvas y residuos)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,6), sharex=True,
                               gridspec_kw={'height_ratios':[3,1]})

for ii in range(len(cpmgs)):
    t = times[ii] * 1000  # convertir a ms
    ydata = cpmgs[ii]/cpmgs[ii][0]
    if substract_offset:
        ydata = ydata - np.mean(ydata[-100:])  # restar offset
    # Valores iniciales para ajuste
    A1_guess = 0.5
    A2_guess = 0.5
    T1_guess = t[int(len(t)*0.3)]
    T2_guess = t[int(len(t)*0.7)]
    y0_guess = 0.0
    p0 = [A1_guess, T1_guess, A2_guess, T2_guess, y0_guess]
    
    # Ajuste biexponencial con offset
    popt, _ = curve_fit(biexp_offset, t, ydata, p0=p0, bounds=(0, np.inf))
    A1, T1, A2, T2, y0 = popt
    
    # Ordenar T1short y T1long
    if T1 < T2:
        T1short, T1long = T1, T2
        P1short = A1 / (A1 + A2)
        P2 = A2 / (A1 + A2)
    else:
        T1short, T1long = T2, T1
        P1short = A2 / (A1 + A2)
        P2 = A1 / (A1 + A2)
    
    # Guardar resultados
    T1short_list.append(T1short)
    T1long_list.append(T1long)
    P1short_list.append(P1short)
    offset_list.append(y0)
    
    # Color para la curva
    color = colors[ii % len(colors)]
    
    # Graficar señal con color Set1 + alpha
    ax1.plot(t, ydata, 'o-', color=color, alpha=0.3, label=f"{samples[ii]}")
    # Graficar ajuste en negro
    yfit = biexp_offset(t, *popt)
    ax1.plot(t, yfit, 'k--', lw=1)
    
    # Graficar residuos
    ax2.plot(t, ydata - yfit, 'o-', color=color, alpha=0.3)
    ax2.axhline(0, color='gray', lw=1, ls='--')
    # Añadir info a la leyenda
    label_text = (f"T2={T1short:.2f} ms ({P1short*100:.0f} %), "
                  f"T2={T1long:.2f} ms ({P2*100:.0f} %)")
    ax1.plot([], [], ' ', label=label_text)  # truco para leyenda solo texto

# Configuración estética
ax1.set_yscale('log')
ax1.set_ylabel("Normalized signal", fontsize=font_size)
ax2.set_xlabel("Time [ms]", fontsize=font_size)
ax2.set_ylabel("Residuals", fontsize=font_size)
ax1.tick_params(axis='both', which='major', labelsize=font_size)
ax2.tick_params(axis='both', which='major', labelsize=font_size)
ax1.legend(fontsize=font_size, loc='upper right')
ax1.set_title("Biexponential fits with offset", fontsize=font_size)
plt.tight_layout()
plt.show()

# %%
