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
expns, sample = [42, 34], 'Li dendrites + KBr'
labels = [r'$\mu$w OFF', r'$\mu$w ON']
npts = 38 # puntos a promediar
nmeans = 5

absolute= False
autoph = False
save = False

path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-12-23_3.2mm_LiMetal/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Overhauser\analysis\2025-12_Enhancement_vs_D1/"
muestra = ""


#=====================================================================
# Ajuste de espectro antes del experimento 1D (before)
#=====================================================================
fig_amps, ax_amps = plt.subplots()
AMPS = []
FIDS = []
colors = ['b', 'r', 'g', 'c', 'm', 'y']
# grafico todos los espectros juntos
# fig_spec, axs_spec = plt.subplots(num=17856, nrows=len(expns), figsize=(6, len(expns)*2))
for jj, expn in enumerate(expns):
    datos = DatosProcesados2D(f'{path}/{expn}/')
    datos.set_fid()
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
    FIDS.append(datos.fid.signal)

    title = f"expn: {expn}\n"
    ax_amps.plot(recycle_delay, np.abs(fid0), 'o-', color=colors[jj], label=labels[jj])
    minimo = recycle_delay[np.argmin(np.abs(fid0))]
    ax_amps.axvline(minimo, color=colors[jj], lw=0.5, ls='--', label=f"RD: {minimo*1e3:.1f} ms")
#%%=====================================================================
nplots = int(npts)
fig_fids, axs_fids = plt.subplots(nrows=nplots, ncols=2, figsize=(6, nplots*2))
for nn in range(len(FIDS)):
    for ii in range(nplots):
        ax_fids = axs_fids[ii, nn]
        label = 'uw OFF' if nn==0 else 'uw ON'
        ax_fids.plot(t*1000, np.abs(FIDS[nn][ii,:]), '-', label=f'{label} \n RD: {recycle_delay[ii]*1000:.3f} ms')
        # ax_fids.plot(t*1000, np.real(FIDS[nn][ii,:]), '-', label=f'{recycle_delay[ii]*1000:.3f} ms')
        # ax_fids.plot(t*1000, np.imag(FIDS[nn][ii,:]), '-')
        ax_fids.set_xlabel("Acquisition time [ms]")
        ax_fids.set_ylabel("abs(FID)[a.u]")
        ax_fids.legend(fontsize=8)
        ax_fids.label_outer()
#%%=====================================================================
# nplots = int(npts)
# fig_fids, axs_fids = plt.subplots(nrows=nplots, ncols=1, figsize=(6, nplots*2))
# for nn in [0,1]:
#     for ii in range(nplots):
#         ax_fids = axs_fids[ii]
#         label = 'uw OFF' if nn==0 else 'uw ON'
#         ax_fids.plot(t*1000, np.abs(FIDS[nn][ii,:])/np.max(np.abs(FIDS[nn][ii,:])), '-', label=f'{label} - RD: {recycle_delay[ii]*1000:.3f} ms')
#         ax_fids.set_xlabel("Acquisition time [ms]")
#         ax_fids.set_ylabel("abs(FID)[norm]")
#         ax_fids.legend()
#         ax_fids.label_outer()
#%%

#=====================================================================
ax = ax_amps
# ax.set_title(r"Signal amplitude \textit{vs} Recycle Delay ")
ax.set_xlabel("Recycle Delay [s]")
ax.set_xscale('log')
ax.legend()
ax.set_ylabel("Signal amplitude [a.u]")
ax.set_yscale('log')

fig, ax = plt.subplots()
enhancement = np.abs(AMPS[1])/np.abs(AMPS[0])
ax.plot(recycle_delay, enhancement, 'ko-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1, None)
ax.grid(True)
ax.set_xlabel("Recycle Delay [s]")
#%%
fig, ax = plt.subplots()
for ii in range(2):
    angle = np.angle(AMPS[ii]) * 180/np.pi
    ax.plot(recycle_delay, angle, 'o-', color=colors[ii], label=labels[ii])
ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel("Recycle Delay [s]")
ax.set_ylabel(r"Phase [$^{\circ}$]")
ax.legend()

# %%

fig, axs = plt.subplots(nrows=2, ncols=1)
ax = axs[0]
norm = np.abs(AMPS[0][-1])
for jj in range(2):
    amp = np.abs(AMPS[jj])
    ax.plot(recycle_delay, amp/norm, 'o-', color=colors[jj], label=labels[jj])
ax.set_xlabel("Recycle Delay [s]")
ax.set_ylabel("Signal amplitude [a.u]")
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.legend()
ax = axs[1]
enhancement = np.abs(AMPS[1])/np.abs(AMPS[0])
ax.plot(recycle_delay, enhancement, 'ko-')
ax.set_xlabel("Recycle Delay [s]")
ax.set_ylabel("Signal ON/OFF")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1, None)
ax.grid(True)
ax.legend()


# %%
