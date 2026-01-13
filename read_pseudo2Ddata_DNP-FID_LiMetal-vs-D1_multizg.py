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




# metal
# expns = np.concatenate([np.arange(1, 60), np.arange(61, 100)])  # directorio de datos
expns, sample = np.arange(20,29), 'Li dendrites + KBr'
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phase90then0/"
Ns = np.arange(2,11)
labels = [f'N = {N}' for N in Ns]

expns, sample = np.arange(30,39), 'Li dendrites + KBr'
expns[0] = 20
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phase90_then-loop-of-0-180/"
Ns = np.arange(2,11)
labels = [f'N = {N}' for N in Ns]

expns, sample = np.arange(40,54), 'Li dendrites + KBr'
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/phases0_vsP1/"
Ns = np.ones_like(expns) * 2

npts = 38 # puntos a promediar
nmeans = 20

absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-12-30_3.2mm_LiMetal-NO-DNP/"
# directorio de guradado
muestra = ""


#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
AMPS = []
p1s = []
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.set_fid()
    p1s.append(datos.acqus.P1)
    # extraigo los primeros npts de la fid
    fid0 = np.mean(datos.fid.real[:,:nmeans], axis=1) + 1j*np.mean(datos.fid.imag[:, :nmeans], axis=1)
    fid0 = fid0[:npts]

    if jj==1:
        slice = 30
        fig0, ax0 = plt.subplots()
        t = datos.acqus.DW * np.arange(datos.fid.real[0,:].size)
        ax0.plot(t*1000, datos.fid.real[slice,:], 'o-', label='Real')
        ax0.plot(t*1000, datos.fid.imag[slice,:], 'o-', label='Imaginary')
        ax0.axvspan(t[0]*1000, t[nmeans]*1000, color='gray', alpha=0.5, label="averaged")
        # ax0.set_xlim(-t[1]*1000, t[4*nmeans]*1000)
        ax0.set_xlabel("Acquisition time [ms]")
        ax0.legend()

    datos.vdlist_path = f"lists/vd/sm2974.D1"
    recycle_delay = datos.get_vdlist(units='s')
    recycle_delay = recycle_delay[:npts]
    AMPS.append(fid0)
    

    title = f"expn: {expn}\n"
    #ax_amps.plot(recycle_delay, np.abs(fid0), 'o-', label=labels[jj])
    
#%%
fig_amps, ax_amps = plt.subplots()
norm = np.abs(AMPS[0][-1])
ax = ax_amps
labels = [f'P1 = {p1}' for p1 in p1s]
for ii in range(Ns.size):
    fid0 = AMPS[ii]
    ax.plot(recycle_delay, np.abs(fid0)/norm, 'o-', label=labels[ii])
#minimo = recycle_delay[np.argmin(np.abs(fid0))]
# ax_amps.axvline(minimo, lw=0.5, ls='--', label=f"RD: {minimo*1e3:.1f} ms")    
# ax.set_title(r"Signal amplitude \textit{vs} Recycle Delay ")
ax.set_xlabel("Recycle Delay [s]")
ax.set_xscale('log')
ax.legend()
ax.set_ylabel("Signal amplitude [a.u]")
ax.set_yscale('log')
ax.grid(True)


fig, ax = plt.subplots()
for ii in range(Ns.size):
    angle = np.angle(AMPS[ii]) * 180/np.pi
    ax.plot(recycle_delay, angle, 'o-', label=labels[ii])
ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel("Recycle Delay [s]")
ax.set_ylabel(r"Phase [$^{\circ}$]")
ax.legend()

for ii in range(expns.size):
    fid0 = AMPS[ii]
    N = Ns[ii]
    p1 = p1s[ii]
    header = "recycle delay [s]\tABS\tphase\tReal\tImag"
    data = np.array([recycle_delay, np.abs(fid0)/norm, np.angle(fid0), np.real(fid0), np.imag(fid0)])
    np.savetxt(f'{savepath}/{N}pulses_p1{p1:.4f}.dat', data, header=header)
np.savetxt(f'{savepath}/P1_list.dat', p1s, header="P1 list in us")

#### ESTA ES PA HACER ENHANCEMENTS
# fig, ax = plt.subplots()
# enhancement = np.abs(AMPS[1])/np.abs(AMPS[0])
# ax.plot(recycle_delay, enhancement, 'ko-')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(1, None)
# ax.grid(True)
# ax.set_xlabel("Recycle Delay [s]")

# %%

# fig, axs = plt.subplots(nrows=2, ncols=1)
# ax = axs[0]
# norm = np.abs(AMPS[0][-1])
# for jj in range(2):
#     amp = np.abs(AMPS[jj])
#     ax.plot(recycle_delay, amp/norm, 'o-', color=colors[jj], label=labels[jj])
# ax.set_xlabel("Recycle Delay [s]")
# ax.set_ylabel("Signal amplitude [a.u]")
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.grid(True)
# ax.legend()
# ax = axs[1]
# enhancement = np.abs(AMPS[1])/np.abs(AMPS[0])
# ax.plot(recycle_delay, enhancement, 'ko-')
# ax.set_xlabel("Recycle Delay [s]")
# ax.set_ylabel("Signal ON/OFF")
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim(9, None)
# ax.grid(True)


