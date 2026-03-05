# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022


@author: Santi



CORREGIR: TIEMPOS Y FASES
"""

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# directorio de datos
expns = np.arange(50, 220, 20)
npts = 50

absolute= False
autoph = False
substract_baseline = False
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-10-01_ccATMC_Rui-R3_NMC-Cu_PC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-10_R3/"
muestra = "7Li_cellR3-PCprotocol"

save = False
plotRange = [400, 125]
# rango de integracion
ppmRange = plotRange# [290, 230]

# #### diamagnetic
# plotRange = [100, -100]
# # rango de integracion
# ppmRange = [50,-50]

# tiempo inicial (para restar y tener tiempos relativos)
datos = DatosProcesados(f'{path}/10/')
t_0 = datos.acqus.dic['DATE_START']

colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))
fig_nut, axs_nut = plt.subplots(num=178156, nrows=2, figsize=(6, 4))

spec_list = []
signal_list = []

signals = np.zeros(expns.size)
tau = np.zeros(expns.size)
ppm_of_max = np.zeros(expns.size)
ppm_mean = np.zeros(expns.size)
phases = np.zeros(expns.size)
for jj, expn in enumerate(expns):
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    # rango de integracion
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    ppmAxis = datos.espectro.ppmAxis

    spec_time = datos.acqus.dic['DATE_START']
    tau[jj] = (spec_time - t_0)/3600  # convert to hours

    
    signal_re = datos.Integrar2D(ppmRange=ppmRange, imaginary=False)[:npts]
    signal_im = datos.Integrar2D(ppmRange=ppmRange, imaginary=True)[:npts]
    tp = np.arange(npts)*2 # us


    axs_nut[0].plot(tp, signal_re/np.max(signal_re), 'o-', label=f'nexp: {expn}')
    axs_nut[1].plot(tp, signal_im, 'o-', label=f'nexp: {expn}')

    ax_spec.plot(datos.espectro.ppmAxis, datos.espectro.real[2, :], label=f'nexp: {expn}')

    spec = np.array([datos.espectro.ppmAxis, datos.espectro.real[2, :], datos.espectro.imag[2, :]]).T
    spec_list.append(spec)
    signal_list.append(signal_re + 1j*signal_im)

#%%
# ============================================================
# Modelo: damped sinewave
# ============================================================
def damped_sine(t, A, T, f, phi, C):
    return A * np.exp(-t / T) * np.sin(2 * np.pi * f * t + phi) + C


# ============================================================
# Leer archivo externo (tiempos + amplitudes)
# ============================================================
path = r"c:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-10_R3\1dspectra\info.txt"

data_ext = np.loadtxt(path)
tau_ext  = data_ext[:, 1]
amp_ext  = data_ext[:, 2]


# ============================================================
# Arrays para guardar parámetros
# ============================================================
A_list   = []
f_list   = []
phi_list = []
T_list   = []
R2_list  = []


# ============================================================
# Loop principal: ajuste + gráficos individuales
# ============================================================
for i, (signal, spec) in enumerate(zip(signal_list, spec_list)):

    y = np.real(signal)

    # --- guesses iniciales ---
    A0 = (y.max() - y.min()) / 2
    T0 = (tp.max() - tp.min()) / 2
    f0 = 1.0 / (tp[1] - tp[0]) / 10
    phi0 = 0.0
    C0 = y.mean()

    p0 = [A0, T0, f0, phi0, C0]

    # --- ajuste ---
    try:
        popt, pcov = curve_fit(damped_sine, tp, y, p0=p0)
        y_fit = damped_sine(tp, *popt)

        A, T, f, phi, C = popt

        # --- goodness of fit: R^2 ---
        ss_res = np.sum((y - y_fit)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        R2 = 1.0 - ss_res / ss_tot

    except RuntimeError:
        A, T, f, phi = np.nan, np.nan, np.nan, np.nan
        R2 = np.nan
        y_fit = None

    # --- guardar parámetros ---
    A_list.append(A)
    T_list.append(T)
    f_list.append(f)
    phi_list.append(phi)
    R2_list.append(R2)

    # ========================================================
    # Figura señal vs tp
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(tp, y, 'o', ms=4, label='Real(signal)')
    ax.plot(tp, np.imag(signal), 'o', ms=4, label='iMAG(signal)')
    if y_fit is not None:
        ax.plot(tp, y_fit, '-', lw=2, label='Damped sine fit')

    ax.set_xlabel('tp')
    ax.set_ylabel('Signal (real)')
    ax.set_title(f'Signal {i}   (R² = {R2:.3f})')
    ax.grid(True)
    ax.legend()

    # ========================================================
    # Inset: espectro NMR
    # ========================================================
    ppm = spec[:, 0]
    spec_real = spec[:, 1]
    spec_imag = spec[:, 2]

    ax_in = inset_axes(ax, width="40%", height="40%", loc='upper right')
    ax_in.plot(ppm, spec_real, lw=1)
    ax_in.plot(ppm, spec_imag, lw=1)
    ax_in.set_xlim(ppm.max(), ppm.min())  # ppm decreciente

    ax_in.set_xlabel('ppm', fontsize=8)
    ax_in.set_ylabel('Intensity', fontsize=8)
    ax_in.tick_params(axis='both', labelsize=8)

    plt.tight_layout()
    plt.show()


# ============================================================
# Conversión a arrays
# ============================================================
A_list   = np.array(A_list)
T_list   = np.array(T_list)
f_list   = np.array(f_list)
phi_list = np.array(phi_list)
R2_list  = np.array(R2_list)


# ============================================================
# Subplots finales: parámetros vs tau
# ============================================================
fig, axs = plt.subplots(
    nrows=6,
    ncols=1,
    figsize=(6, 10),
    sharex=True
)

# --- Datos externos ---
axs[0].plot(tau_ext, amp_ext, 'o-')
axs[0].set_ylabel('External amplitude')
axs[0].set_title('External data')
axs[0].grid(True)

# --- Amplitud ---
axs[1].plot(tau, A_list, 'o-')
axs[1].set_ylabel('Amplitude')
axs[1].grid(True)

# --- Frecuencia ---
axs[2].plot(tau, f_list, 'o-')
axs[2].set_ylabel('Frequency')
axs[2].grid(True)

# --- Fase ---
axs[3].plot(tau, phi_list, 'o-')
axs[3].set_ylabel('Phase (rad)')
axs[3].grid(True)

# --- Tiempo característico del damping ---
axs[4].plot(tau, T_list, 'o-')
axs[4].set_ylabel('T (damping)')
axs[4].grid(True)

# --- Goodness of fit ---
axs[5].plot(tau, R2_list, 'o-')
axs[5].set_ylabel('R²')
axs[5].set_xlabel('tau')
axs[5].grid(True)

plt.tight_layout()
plt.show()

# %%
# # --- Amplitud ---
# fig, ax = plt.subplots(figsize=(4, 6))
# ax.plot(A_list, f_list, 'o')
# ax.set_xlabel('Amplitude')
# ax.set_ylabel('Frequency')
# ax.grid(True)

#%%
for n in range(len(spec_list)):
    spec = spec_list[n]

    ppm = spec[:, 0]
    spec_real = spec[:, 1]
    spec_imag = spec[:, 2]
    spec_complex = spec_real + 1j*spec_imag

    spec90 = spec_complex * np.exp(-1j*np.pi/2)

    spec90_real = np.real(spec90)
    spec90_imag = np.imag(spec90)

    plt.figure(n)
    plt.plot(ppm, spec_real, label='real', lw=3)
    plt.plot(ppm, spec_imag, label='imag', lw=3)
    plt.plot(ppm, spec90_real, label='real 90deg', lw=0.5)
    plt.plot(ppm, spec90_imag, label='imag 90deg', lw=0.5)
    plt.legend()
# %%
