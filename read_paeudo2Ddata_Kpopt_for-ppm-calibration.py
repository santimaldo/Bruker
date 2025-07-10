# -*- coding: utf-8 -*-
"""
Simplificado para encontrar la posición en ppm del máximo de cada espectro,
refinando la posición con interpolación de la derivada del módulo,
con ajuste lineal y gráficos correspondientes, y gráficos individuales
de derivada e interpolación alrededor del máximo.

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Autophase import autophase
import re
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

################## Functions ###########################

def read_bsms_field(path_archivo):
    """
    Lee un archivo con formato específico y devuelve el array 'x' como lista de floats.
    """
    with open(path_archivo + 'Klog_opt', 'r') as f:
        for linea in f:
            if linea.startswith("x="):
                contenido = linea.strip().split("=", 1)[1]
                x = eval(contenido)
                return np.array(x)
    raise ValueError("No se encontró una línea que empiece con 'x='.")


def calc_derivative_and_interp(ppm_axis, spectrum_mod, window_pts=3, interp_points=500):
    """
    Calcula la derivada y su interpolación en la ventana alrededor del máximo.

    Retorna:
    - ppm_window: eje ppm en ventana
    - deriv: derivada discreta en ventana
    - ppm_interp: eje ppm interpolado
    - deriv_interp: derivada interpolada
    - idx_max: índice del máximo bruto dentro del espectro completo
    """
    idx_max = np.argmax(spectrum_mod)

    # Definir ventana alrededor del máximo
    start = max(idx_max - window_pts, 0)
    end = min(idx_max + window_pts + 1, len(spectrum_mod))

    ppm_window = ppm_axis[start:end]
    mod_window = spectrum_mod[start:end]

    # Calcular derivada discreta (dmod/dppm)
    deriv = np.gradient(mod_window, ppm_window)

    # Interpolación spline cúbica para la derivada
    interp_func = interp1d(ppm_window, deriv, kind='cubic')
    ppm_interp = np.linspace(ppm_window.min(), ppm_window.max(), interp_points*ppm_axis.size)
    deriv_interp = interp_func(ppm_interp)

    return ppm_window, deriv, ppm_interp, deriv_interp, idx_max


def refine_peak_position_from_deriv(ppm_interp, deriv_interp, ppm_max):
    """
    Encuentra el punto de cruce por cero más cercano al máximo original
    a partir de la derivada interpolada.

    Retorna:
    - ppm_refined: posición refinada en ppm
    """
    sign_changes = np.where(np.diff(np.sign(deriv_interp)))[0]

    ppm_refined = 0.5*(ppm_interp[sign_changes]+ppm_interp[sign_changes+1]) 
    return float(ppm_refined)

################## end Functions #######################

# Parámetros
expns = [45]
absolute = False
autoph = False
path  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-06-16_3.2mm_IMECdendrites/"
plotRange = [40, -40]

# Gráficos
fig_spec, ax_spec = plt.subplots(num=17856)
fig_ppm, ax_ppm = plt.subplots(num=29912)

# Para gráficos de derivadas/interpolaciones:
deriv_fig, deriv_axes = None, []

for expn in expns:
    ppm_max_positions = []
    ppm_refined_positions = []
    deriv_data = []  # para guardar datos de derivadas e interpolación para graficar luego

    path_2D = f"{path}/{expn}/"
    datos = DatosProcesados2D(path_2D, read_pp=False)
    datos.espectro.ppmSelect(plotRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    speci = datos.espectro.imag

    # Procesamiento espectros
    for kk in range(spec.shape[0]):
        spec1d = spec[kk, :]
        speci1d = speci[kk, :]

        if absolute:
            abs_spec = np.abs(spec1d + 1j * speci1d)
            spec_mod = abs_spec
        elif autoph:
            spec1d = ng.proc_autophase.autops(spec1d + 1j * speci1d, "acme").real
            spec1d = ng.process.proc_bl.cbf(spec1d, last=100)
            spec_mod = np.abs(spec1d)
        else:
            spec_mod = spec1d 

        ppm_max = ppmAxis[np.argmax(spec_mod)]
        ppm_max_positions.append(ppm_max)

        # Obtener datos de derivada e interpolación
        ppm_window, deriv, ppm_interp, deriv_interp, idx_max = calc_derivative_and_interp(ppmAxis, spec_mod)
        ppm_refined = refine_peak_position_from_deriv(ppm_interp, deriv_interp, ppm_max)
        ppm_refined_positions.append(ppm_refined)

        deriv_data.append((ppm_window, deriv, ppm_interp, deriv_interp))

        ax_spec.plot(ppmAxis, spec_mod, label=f'Spectrum {kk}' if kk==0 else "", alpha=0.7)

    bsms_field = read_bsms_field(path_2D)
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax_spec.axhline(0, color='k')
    ax_spec.set_title('Espectros módulo')
    ax_spec.legend()

    # Mostrar resultados
    print("bsms_field | ppm máximo (bruto) | ppm máximo (refinado)")
    for b, p_max, p_ref in zip(bsms_field, ppm_max_positions, ppm_refined_positions):
        print(f"{b:.4f}     {p_max:.4f}     {p_ref:.4f}")

    # Gráfico ppm máximo vs bsms_field (refinado)
    ax_ppm.plot(bsms_field, ppm_refined_positions, 'o', color='darkred', label='ppm refinado')
    ax_ppm.set_xlabel("BSMS field (arb. units)")
    ax_ppm.set_ylabel("ppm máximo (refinado)")
    ax_ppm.set_title("ppm del máximo refinado vs campo BSMS")
    ax_ppm.grid(True)

    # --- Ajuste lineal: ppm = a * bsms + b
    def linear(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(linear, bsms_field, ppm_refined_positions)
    a_fit, b_fit = popt
    ppm_fit = linear(bsms_field, *popt)
    residuals = np.array(ppm_refined_positions) - ppm_fit

    print(f"\nAjuste lineal: ppm = {a_fit:.6f} × bsms + {b_fit:.6f}")


    # --- Gráficos de derivada/interpolación para cada espectro ---
    n = len(deriv_data)
    deriv_fig, deriv_axes = plt.subplots(n, 1, figsize=(6, 3*n), sharex=False, num=556677)
    if n == 1:
        deriv_axes = [deriv_axes]

    for i, (ax, (ppm_w, deriv_w, ppm_i, deriv_i)) in enumerate(zip(deriv_axes, deriv_data)):
        ax.plot(ppm_i, deriv_i, 'o', label='Derivada interpolada', markersize=1)
        ax.plot(ppm_w, deriv_w, 'o', label='Derivada discreta', markersize=8)
        ax.axhline(0, color='k', linestyle='--')
        ax.set_xlim(ppm_w.min(), ppm_w.max())
        ax.set_ylabel(f'bsms_field {bsms_field[i]}')
        ax.legend()
        ax.grid(True)

    deriv_axes[-1].set_xlabel('ppm')



    # --- Figura con ajuste y residuales
    fig_fit, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6),
                                       gridspec_kw={"height_ratios": [3, 1]}, num=122334)
    ax1.plot(bsms_field, ppm_refined_positions, 'o', label='Datos', color='darkred')
    ax1.plot(bsms_field, ppm_fit, '-', label=f'Ajuste: y = {a_fit:.3e}·x + {b_fit:.3f}', color='black')
    ax1.set_ylabel('ppm máximo refinado')
    ax1.legend()
    ax1.grid(True)

    ax2.axhline(0, color='gray', linestyle='--')
    ax2.plot(bsms_field, residuals, 'o', color='gray')
    ax2.set_xlabel('BSMS field (arb. units)')
    ax2.set_ylabel('Residuos')
    ax2.grid(True)

plt.tight_layout()
plt.show()

#%%
def get_SR(bsms, reference_ppm):
    """
    Calcula el offset para el campo BSMS.
    """
    # bsms=-2400
    # reference_ppm = 1.89
    delta0 = reference_ppm
    delta = linear(bsms, *popt)
    sf = (datos.procs.dic["SF"]) * 1e6
    bf1 = datos.acqus.dic["BF1"] * 1e6
    sr = sf - bf1 
    sfo1 = datos.acqus.SFO1
    sf0 =  sf + sfo1 *(delta - delta0)

    sr0 = sf0 - bf1 
    print("delta0: ", delta0)
    print("delta: ", delta)
    print("SFO1: ", sfo1)
    print("SF: ", sf)
    print("BF1: ", bf1)
    print("SR: ", sr)
    print("SR0: ", sr0)
    print("bsms: ", bsms)
    print("ppm referencia: ", reference_ppm)



    print("era: ", sr)
    print("ahora: ", sr0)
    print("se corrio en ppm: ", (sr - sr0)/sfo1)
    return sr0

#%%
srs = []
for bsms in bsms_field:
    get_SR(bsms, 1.89)  # Ajusta el valor de referencia_ppm según sea necesario
    srs.append(get_SR(bsms, 1.89))
plt.figure(num=123456)
plt.plot(bsms_field, srs, 'o', color='darkblue', label='SR calculado')

# %%
