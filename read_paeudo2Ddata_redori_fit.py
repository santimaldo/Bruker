# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 2025

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
import scipy.integrate as integrate
import re
from redor import *
from VoigtFit import VoigtFit

#################### functions ################################################


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def label_curve(ax, x, y, label, idx, offset=(0, 0), **kwargs):
    """Anota una curva en el punto `idx`, sin rotación."""
    x0 = x[idx] + offset[0]
    y0 = y[idx] + offset[1]

    ax.text(x0, y0, label,
            rotation=0,
            ha='center', va='center',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.5),
            **kwargs)

############################################################

# Directorio de datos
expn = 44
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-06-21_PEO-solid-electrolyte/"
# Directorio de guardado
savepath = r"C:/"
muestra = ""
save = False
plot_individual_pairs = False  # Activar/desactivar gráficos por par

plotRange = [-52, -70]
peaks = [-40, -80]
range_of_peaks_to_save = [-50, -65]  # Rango de ppm para guardar los picos

last_echo_number = 25

#=====================================================================
# Lectura del experimento 2D
#=====================================================================

path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(path_2D, read_pp=False)
datos.espectro.ppmSelect(plotRange)
ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real


m1_amplitude = 227.150490 # initial guess for m1 amplitude
m2_amplitude = 100.065480
# Graficar los espectros 1D
Signals = np.array([])
Signals_err = np.array([])

for kk in range(2*last_echo_number):
    ydata = spec[kk, :]
    xdata = ppmAxis
    vfit=VoigtFit(xdata,
              ydata,
              Npicos=2,
              ajustar=True,
              amplitude=[m1_amplitude, m2_amplitude],
              center=[-59.7286635, 	-59.1561544],
              sigma=[1.15184146, 9.3274e-06],
              gamma=[0.38857500,0.62508785],
              fijar=['center', 'sigma', 'gamma'],
              )
    fig = vfit.plot_ajuste()
    fig.gca().set_title(f"Slice Num. {kk}")

    m1_amplitude = vfit.params["m1_amplitude"].value
    m2_amplitude = vfit.params["m2_amplitude"].value
    signal = m1_amplitude + m2_amplitude

    m1_amplitude_stderr = vfit.params["m1_amplitude"].stderr
    m2_amplitude_stderr = vfit.params["m2_amplitude"].stderr
    signal_stderr = np.sqrt(m1_amplitude_stderr**2 + m2_amplitude_stderr**2)

    Signals = np.append(Signals, signal)
    Signals_err = np.append(Signals_err, 1.96*signal_stderr)
    # ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    # ax_spec.axhline(0, color='k')
    # ax_spec.set_title(f"Echo {kk}")




#=====================================================================
# Calculo y grafico de (S - S0)/S0
#=====================================================================

# Separar S y S0
S = Signals[0::2]
S0 = Signals[1::2]
S_err = Signals_err[0::2]
S0_err = Signals_err[1::2]
N = np.arange(1, len(S) + 1) # number of rotor cycles

# Calcular la razón (S - S0)/S0
spin_speed = 14000  # spinning speed in Hz
recopl_time = N / spin_speed  # recoupling time in seconds

# Graficar S y S0 en bruto
fig_t2, ax_t2 = plt.subplots()
ax_t2.plot(S, 'o-', label='S')
ax_t2.plot(S0, 'o-', label='S0')
ax_t2.set_xlabel("Índice (simula tiempo)")
ax_t2.set_ylabel("S")
ax_t2.grid(True)
ax_t2.legend()

#%% Graficar S y S0 en bruto
fig_redor, ax_redor = plt.subplots()
for internuc_distance in np.arange(0.55, 0.8, 0.05):
    Nsim = np.arange(1, 1.5*N[-1])  # rotor cycles from 1 to 64
    NTr, S_S0 = DeltaS_quadrupolar(internuc_distance, spin_speed, Nsim)
    xvals = NTr * 1000
    ax_redor.plot(xvals, S_S0, 'k')
    label = fr"{internuc_distance*10:.1f} $\AA$"
    label_curve(ax_redor, xvals, S_S0, label, idx=int(1.1*N[-1]), offset=(0, 0.05), fontsize=8)

Necos = 25# for the plot]
S_S0_Necos = S[:Necos]/S0[:Necos]
S_S0_err = S_S0_Necos * np.sqrt( (S_err[:Necos]/S[:Necos])**2 + (S0_err[:Necos]/S0[:Necos])**2 )

ax_redor.errorbar(recopl_time[:Necos]*1000, S_S0_Necos ,
                  yerr=S_S0_err,
                  fmt='o', alpha=0.8, capsize=3)
ax_redor.set_xlabel("Dephasing time (ms)")
ax_redor.set_ylabel("S/S0")
ax_redor.set_ylim(0, 1.1)
ax_redor.grid(True)
ax_redor.legend()


#%% Plot S and S0 vs recoupling time
fig_t2, ax_t2 = plt.subplots()
ax_t2.plot(recopl_time[:Necos]*1000, S[:Necos]/S[0],'o-', label='S')
ax_t2.plot(recopl_time[:Necos]*1000, S0[:Necos]/S[0], 'o-', label='S0')
ax_t2.set_xlabel("Dephasing time (ms)")
ax_t2.set_ylabel("Signal (a.u.)")
ax_t2.set_title("S and S0 vs Dephasing time")
ax_t2.grid(True)
ax_t2.legend()
ax_t2.set_xlim(ax_redor.get_xlim())
ax_t2.set_yscale('log')
# ax_t2.set_ylim(0.0, 0.05)
# %%
