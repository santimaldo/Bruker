# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:09:33 2022

@author: Santi

Extrae multiples espectros adquiridos en Bruker Avance II
"""

# import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from Datos import *
from Espectro import autophase
from scipy.stats import linregress
import matplotlib.ticker as ticker


# Nmuestra es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
Nmuestra = 10
# directorio de datos:
fecha = '10-11'  # 'MM-DD'
expnums = np.arange(20, 30)
# rango de integracion
ppmRange = [-1, 1]
save = True

min_gp = 5


#-------------------- directorios
path_local = "S:/CNEA/Glicerol-Agua/116MHz"
path_bruker = f"/2022-{fecha}_Diff_Silica_Agua-Glicerol-LiCl/"
path = path_local + path_bruker
# directorio de guradado
savepath = "S:/CNEA/Glicerol-Agua/analisis/"


# -----------------------------------------------
# datos de la muestra
muestra = f"M{Nmuestra}"
N = Nmuestra-10
# matriz es el medio: Bulk, Q30 o Q3
if N < 3:
    matriz = 'Q3'
elif N < 6:
    matriz = 'Q30'
elif N < 9:
    matriz = 'Bulk'
# pc es el porcentaje de glicerol: 50%, 70% o 90%
if N % 3 == 0:
    pc = 50
elif N % 3 == 1:
    pc = 70
elif N % 3 == 2:
    pc = 90
msg = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}"
print(msg)

# %%
gplist = []
intensidades = []
intensidadesFID = []
fig1d, axs = plt.subplots(nrows=1, ncols=2)  # create figure & 1 axis
for nn in range(len(expnums)):
    expn = expnums[nn]
    print(expn)
    # extraigo:
    datos = DatosProcesados(f'{path}/{expn}/')
    nucleo = datos.nucleo
    gplist.append(datos.acqus.gp)
    datos.espectro.ppmSelect(ppmRange)
    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis

    integral = scipy.integrate.simps(re, x=-ppmAxis)
    intensidades.append(integral)

    # calculo FID
    datos.set_fid()
    timeAxis = datos.fid.timeAxis
    fid = datos.fid.real
    fid = np.abs(datos.fid.real + 1j*datos.fid.imag)
    npts = 20
    intensidadFID = np.abs(np.sum(fid[0:npts]))
    intensidadesFID.append(intensidadFID)

    # guardo:
    # header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
    # dataexport = np.array([ppmAxis, re, im, re, im]).T
    # filename = f'{savepath}/{nucleo}_{muestra}.dat'
    # # np.savetxt(filename, dataexport, header=header)

    # grafico para ver:
    print('graficando...', nucleo, muestra)
    axs[0].plot(ppmAxis, re, linewidth=2)
    axs[0].set_xlabel(f"{nucleo} NMR Shift [ppm]")
    axs[0].set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    axs[1].plot(timeAxis*1000, fid, 'o-', linewidth=2)
    axs[1].set_xlabel(f"time [ms]")
    if (datos.acqus.gp == min_gp):
        spec1d = re
        spec1d_im = im

bigDelta = datos.acqus.D[20]  # ms
delta = datos.acqus.P[30]*2 * 1e-6  # ms
titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n" \
         fr"$\Delta = {bigDelta*1e3}$ ms, $\delta = {delta*1e3}$ ms"
fig1d.suptitle(titulo)
# %%
# calculo bvalue:
# gplist = np.loadtxt(f"{path}gp_list.dat")[:,1]
gplist = np.array(gplist)
bigDelta = datos.acqus.D[20]  # s
delta = datos.acqus.P[30]*2 * 1e-6  # s
gamma = 103.962e6  # rad/(s*T) ---> 7Li
g0 = 12  # T/m
bvalue = (gamma*gplist/100*g0*delta)**2*(bigDelta-delta/3) * 1e-9

# Senal:
signal = np.array(intensidades)
# S = np.array(intensidadesFID)

# Ajuste lineal----------------------------
# Reordeno:
signal = signal[bvalue[:].argsort()]
bvalue = bvalue[bvalue[:].argsort()]
# defino variables de ajuste
y = np.log(signal)
x = bvalue
slope, intercept, r, p, se = linregress(x, y)
yfit = slope*x+intercept
signal_fit = np.exp(yfit)
residuals = signal - signal_fit
# Resultados del ajuste
D = -slope
uD = se
r_squared = r**2
msg = f"Ajuste lineal de Difusion:\n \
  D =  {D:.8f} 10^-9 m^2/s\n \
  Rsquared = {r_squared:.6f}"
print(msg)


# %% Grafico
smax = np.max(signal)
signal = signal/smax
signal_fit = signal_fit/smax
bvalue_fit = bvalue
residuals = residuals/smax

fig, axs = plt.subplots(2, 1, figsize=(6, 7))
# -------------
titulo = f"muestra: M{Nmuestra}, {pc}% glicerol, {matriz}\n" \
         fr"$\Delta = {bigDelta*1e3}$ ms, $\delta = {delta*1e3}$ ms"
axs[0].set_title(titulo)
axs[0].plot(bvalue, signal, 'ko')
axs[0].plot(bvalue_fit, signal_fit, 'r-')
text = f"$D =$ {D:.4f} $10^{{-9}}m^2/s$ \n " \
       f"$u_D =$ {uD:.1g} $10^{{-9}}m^2/s$ \n " \
       f"$r^2 =$ {r_squared:.5g}"
axs[0].text(bvalue[-1]*0.6, 0.7*signal[0], text,
            multialignment="left")
axs[0].set(xlabel=r'$b_{value} [10^9 s/m^2]$', ylabel=r'$S$')
axs[0].set_yscale('log')
axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1g'))
axs[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1g'))
# -------------
axs[1].plot(bvalue, residuals, 'ko')
axs[1].axhline(0, color='k', linestyle='--')
axs[1].set(xlabel=r'$b_{value} [10^9 s/m^2]$', ylabel=r'Residuos')
axs[1].set_yscale('symlog')
ylim = 10**np.ceil(np.log10(np.max(abs(residuals))))
axs[1].set_ylim([-ylim, ylim])

# %%
# guardo data:
if save:
    filename0 = f"{muestra}_{pc}pc_{matriz}"

    filename = f'{savepath}/datos_Diff/figuras/{filename0}_Diff.png'
    fig.savefig(filename)   # save the figure to file

    filename = f'{savepath}/datos_Diff/figuras/'\
               f'{filename0}_Diff-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n" \
             f"Rango de integracion: {ppmRange} ppm\n" \
             f"bvalue (10^-9 m^2/s)\t S (norm)"
    Diffdata = np.array([bvalue, signal]).T
    np.savetxt(f"{savepath}/datos_Diff/{filename0}_Diff.dat",
               Diffdata, header=header)

    data = np.array([ppmAxis, spec1d, spec1d_im]).T
    filename = f"{savepath}/datos_Diff/primerEspectro/" \
               f"{filename0}_primerEspectro.dat"
    np.savetxt(filename, data)
