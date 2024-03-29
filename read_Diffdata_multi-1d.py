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


# Nmuestra:  es el n de Mn, ejemplo: M16 ---> Nmuestra = 16
# fecha: correspondiente al directorio de datos.  'MM-DD'
# expn: numero de archivo
# ppmRange: rango de integracion
# bmax: maximo valor de bvalue utilizado para ajustar
#
# info = [muestra, expni, expnf, ppmRange]

# M4 HF LiTFSI
path_bruker = f"2023-10-10_CarbonesCNEA_M4_Diff/"
info = ['7Li_M4-HF-LiTFSI', 10, 25, [-2,2]]




min_gp = 5

gpshape = 'sin'
factor_b = 1

FIDsignal = True
NptsFid = 
ABSsignal = False  # abs del espectro
centrado_en_maximo = False

save = True
save1d = False
muestra, expni, expnf, ppmRange = info
expnums = np.arange(expni, expnf+1)
# expnums = [28]

#-------------------- directorios
# carbones CNEA
path_local = "S:/NMRdata/2018_Carbones_CNEA/"
path = path_local + path_bruker
# directorio de guradado
# savepath_local = "S:/"  # Oficina
# savepath_local = "G:/Otros ordenadores/Oficina/" # Acer
# savepath = f"{savepath_local}Posdoc/CNEA/Glicerol-Agua/analisis/7Li/"

savepath = "S:/Posdoc/CNEA/Carbones/analisis/2023-10_Difusion/datos/Diff/"  # Oficina


# %%
gplist = []
intensidades = []
intensidadesFID = []
color = []
fig1d, axs = plt.subplots(nrows=1, ncols=2)  # create figure & 1 axis
for nn in range(len(expnums)):
    expn = expnums[nn]
    print(expn)
    if True:  # antes era "if expn<15:"
        color.append('k')
    else:
        color.append('b')
    # extraigo:
    datos = DatosProcesados(f'{path}/{expn}/')
    nucleo = datos.nucleo
    gplist.append(datos.acqus.gp)
    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis

    # guardo espectros 1d
    if save1d:
        filename0 = f"{muestra}_{pc}pc_{matriz}"
        data = np.array([ppmAxis, re, im]).T
        filename = f"{savepath}/datos_Diff/Espectros_vs_gpz/" \
                   f"{filename0}_gpz{datos.acqus.gp:0>5.2f}.dat"
        np.savetxt(filename, data)

    # recorto para integrar
    datos.espectro.ppmSelect(ppmRange, centrado_en_maximo=centrado_en_maximo)
    re = datos.espectro.real
    im = datos.espectro.imag
    ppmAxis = datos.espectro.ppmAxis

    # plt.figure(expn)
    # plt.axhline(0)
    # plt.plot(ppmAxis, re, 'k')
    # plt.plot(ppmAxis, im, 'r')
    # plt.plot(ppmAxis, np.abs(re+1j*im), 'orange')
    # plt.title(f"GP: {datos.acqus.gp:.2f} %")

    if ABSsignal:
        spec_integrado = np.abs(re+1j*im)
    else:
        spec_integrado = re
    ancho = np.abs(ppmAxis[0]-ppmAxis[-1])
    integral = scipy.integrate.simps(spec_integrado, x=-ppmAxis) / ancho
    intensidades.append(integral)

    # calculo FID
    datos.set_fid()
    timeAxis = datos.fid.timeAxis
    fid = datos.fid.real
    fid = np.abs(datos.fid.real + 1j*datos.fid.imag)
    intensidadFID = np.sum(fid[0:NptsFid])
    intensidadesFID.append(intensidadFID)

    # guardo:
    # header = "ppmAxis\t real (norm)\t imag (norm)\t real \t imag"
    # dataexport = np.array([ppmAxis, re, im, re, im]).T
    # filename = f'{savepath}/{nucleo}_{muestra}.dat'
    # # np.savetxt(filename, dataexport, header=header)

    # grafico para ver:
    print('graficando...', nucleo, muestra)
    axs[0].plot(ppmAxis, spec_integrado, linewidth=2, color=color[nn])
    axs[0].set_xlabel(f"{nucleo} NMR Shift [ppm]")
    # axs[0].set_xlim([np.max(ppmAxis), np.min(ppmAxis)])
    axs[1].plot(timeAxis*1000, fid, 'o-', linewidth=2)
    axs[1].set_xlabel(f"time [ms]")
    if (datos.acqus.gp == min_gp):
        spec1d = re
        spec1d_im = im

bigDelta = datos.acqus.D[20]  # ms
delta = datos.acqus.P[30]*2 * 1e-6  # ms
titulo = f"muestra: {muestra}\n" \
         fr"$\Delta = {bigDelta*1e3}$ ms, $\delta = {delta*1e3}$ ms"
fig1d.suptitle(titulo)
# %%
# calculo bvalue:
# gplist = np.loadtxt(f"{path}gp_list.dat")[:,1]
bigDelta = datos.acqus.D[20]  # s
delta = datos.acqus.P[30]*2 * 1e-6  # s
gamma = 103.962e6  # rad/(s*T) ---> 7Li
g0 = 12  # T/m
gplist = np.array(gplist)*g0/100
if 'sin' in gpshape.lower():
    bvalue = (gamma*gplist*delta/np.pi)**2 * (4*bigDelta-delta) * 1e-9
elif 'rect' in gpshape.lower():
    bvalue = (gamma*gplist*delta)**2 * (bigDelta-delta/3) * 1e-9
bvalue = factor_b * bvalue
# Senal:

if FIDsignal:
    signal = np.array(intensidadesFID)
else:
    signal = np.array(intensidades)
# Ajuste lineal----------------------------
# Reordeno:
signal = signal[bvalue[:].argsort()]
color = np.array(color)[bvalue[:].argsort()]
gplist = gplist[bvalue[:].argsort()]
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
titulo = f"muestra: {muestra}\n" \
         fr"$\Delta = {bigDelta*1e3}$ ms, $\delta = {delta*1e3}$ ms"
axs[0].set_title(titulo)
# axs[0].plot(bvalue, signal, 'ko')
axs[0].scatter(bvalue, signal, c=color)
axs[0].plot(bvalue_fit, signal_fit, 'r-')
text = f"$D =$ {D:.4f} $10^{{-9}}m^2/s$ \n " \
       f"$u_D =$ {uD:.1g} $10^{{-9}}m^2/s$ \n " \
       f"$r^2 =$ {r_squared:.5g}"
axs[0].text(bvalue[-1]*0.6, 0.5*(max(signal)-min(signal))+min(signal), text,
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
    filename0 = f"{muestra}"

    # filename = f'{savepath}/datos_Diff/figuras/{filename0}_Diff.png'
    filename = f'{savepath}/{filename0}_Diff.png'
    fig.savefig(filename)   # save the figure to file

    # filename = f'{savepath}/datos_Diff/figuras/'\
    filename = f'{savepath}/{filename0}_Diff-RegionIntegracion.png'
    fig1d.savefig(filename)   # save the figure to file

    header = f"Archivo de datos: {path_bruker}\n" \
             f"Rango de integracion: {ppmRange} ppm\n" \
             f"bvalue (10^-9 m^2/s)\t S (norm)\n" \
             f"gpz (%)"
    Diffdata = np.array([bvalue, signal, gplist]).T
    # np.savetxt(f"{savepath}/datos_Diff/{filename0}_Diff.dat",
    np.savetxt(f"{savepath}/{filename0}_Diff.dat",
               Diffdata, header=header)

    # data = np.array([ppmAxis, spec1d, spec1d_im]).T
    # filename = f"{savepath}/datos_Diff/primerEspectro/" \
    #            f"{filename0}_primerEspectro.dat"
    # np.savetxt(filename, data)
plt.show()
