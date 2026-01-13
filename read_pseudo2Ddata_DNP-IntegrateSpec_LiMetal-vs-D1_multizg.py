# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:10:32 2022

@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import VoigtFit as vf


#=====================================================================
# Experimentos
#=====================================================================

expns, sample = np.arange(40, 63), 'Li dendrites + KBr'
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phases0_vsP1/"
Ns = np.ones_like(expns) * 2

# numero de puntos (RDs)
npts = 36

# opciones
absolute = False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-12-30_3.2mm_LiMetal-NO-DNP/"
muestra = ""

# parametros de integracion espectral
plotRange = [350, 320]
ppmIntegrationWidth = 100  # ppm


#=====================================================================
# Loop principal
#=====================================================================

AMPS = []
p1s = []

for jj, expn in enumerate(expns):

    datos = DatosProcesados2D(f'{path}/{expn}/')

    # ----------------------------
    # parametros del experimento
    # ----------------------------
    p1s.append(datos.acqus.P1)

    datos.vdlist_path = f"lists/vd/sm2974.D1"
    recycle_delay = datos.get_vdlist(units='s')
    recycle_delay = recycle_delay[:npts]

    # ----------------------------
    # espectro 1D
    # ----------------------------
    datos.espectro.ppmSelect(plotRange)

    ppm = datos.espectro.ppmAxis
    re  = datos.espectro.real
    im  = datos.espectro.imag

    # ----------------------------
    # integral compleja del espectro
    # ----------------------------
    spec_cplx = re + 1j*im
    spec_roi  = spec_cplx

    signal = integrate.simpson(spec_roi, x=ppm, axis=1)
    signal = signal[:npts]

    AMPS.append(signal)

    # ----------------------------
    # diagnostico visual
    # ----------------------------
    if jj == 1:
        fig, ax = plt.subplots()
        ax.plot(ppm, re[jj]/np.max(re[jj]), 'k', lw=2)
        ax.set_xlim(np.max(ppm), np.min(ppm))
        ax.set_xlabel("NMR shift [ppm]")
        ax.set_ylabel("Norm. intensity")
        ax.legend()
        ax.set_title(f"expn {expn}")


#=====================================================================
# Plots: amplitud
#=====================================================================

fig_amps, ax_amps = plt.subplots()

norm = np.abs(AMPS[0][-1])
labels = [f'P1 = {p1:.3f} µs' for p1 in p1s]

for ii in range(len(AMPS)):
    ax_amps.plot(recycle_delay,
                 np.abs(AMPS[ii]) / norm,
                 'o-', label=labels[ii])

ax_amps.set_xlabel("Recycle Delay [s]")
ax_amps.set_ylabel("Signal amplitude [a.u.]")
ax_amps.set_xscale('log')
ax_amps.set_yscale('log')
ax_amps.grid(True)
ax_amps.legend()


#=====================================================================
# Plots: fase
#=====================================================================

fig_phase, ax_phase = plt.subplots()

for ii in range(len(AMPS)):
    phase = np.angle(AMPS[ii]) * 180 / np.pi
    ax_phase.plot(recycle_delay, phase, 'o-', label=labels[ii])

ax_phase.set_xlabel("Recycle Delay [s]")
ax_phase.set_ylabel(r"Phase [$^\circ$]")
ax_phase.set_xscale('log')
ax_phase.grid(True)
ax_phase.legend()

p1s.sort()
#=====================================================================
# Guardado
#=====================================================================

# for ii in range(len(expns)):

#     signal = AMPS[ii]
#     N = Ns[ii]
#     p1 = p1s[ii]

#     header = "recycle delay [s]\tABS\tphase\tReal\tImag"
#     data = np.array([
#         recycle_delay,
#         np.abs(signal) / norm,
#         np.angle(signal),
#         np.real(signal),
#         np.imag(signal)
#     ])

#     np.savetxt(f'{savepath}/{N}pulses_p1{p1:.4f}.dat',
#                data,
#                header=header)

# np.savetxt(f'{savepath}/P1_list.dat',
#            p1s,
#            header="P1 list in us")


#=====================================================================
# Enhancements (opcional)
#=====================================================================

# fig, ax = plt.subplots()
# enhancement = np.abs(AMPS[1]) / np.abs(AMPS[0])
# ax.plot(recycle_delay, enhancement, 'ko-')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1, None)
# ax.grid(True)
# ax.set_xlabel("Recycle Delay [s]")

# %%
