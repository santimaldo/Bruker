# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Autophase import autophase
from Laplace import ILT
from scipy.optimize import curve_fit

save = False

# ============================================================
# DIRECTORIOS
# ============================================================

expns = np.arange(35, 1200, 10)
samples = [""] * len(expns)
path = r"c:\Users\Santi\Documents\NMRdata\300neo\2026-08-28_ATMC1_Rui-R12_NMC-Cu_CC\\"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2026-08_R12_CC\\"

# ============================================================
# LEER DATOS
# ============================================================

cpmgs, times = [], []

for expn, sample in zip(expns, samples):
    datos = Datos(f"{path}{expn}", set_fid=True)
    echotime = 2 * datos.acqus.dic["D"][20]
    re, im = datos.fid.real, datos.fid.imag
    timeAxis = np.arange(1, re.size + 1) * echotime
    fid, angle = autophase(timeAxis, re + 1j * im, method='minIntImag')
    cpmgs.append(fid.real)
    times.append(timeAxis)

#============================================================
# CPMG: SUPERPOSICIÓN VS TIEMPO
# ============================================================

fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.cm.viridis
norm = plt.Normalize(expns.min(), expns.max())

for ii, expn in enumerate(expns):
    timeAxis = times[ii]
    ydata = cpmgs[ii] / cpmgs[ii][0]
    ax.plot(timeAxis * 1000, ydata, color=cmap(norm(expn)), lw=1.5)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax)
cbar.set_label("Experiment number")

ax.set_xlabel("Time [ms]")
ax.set_ylabel("Normalized CPMG amplitude")
ax.set_title("CPMG decays")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# ============================================================
# ILT
# ============================================================

Npts_L = 100
alpha = 1e0

def T2kernel(x, tt):
    return np.exp(-x / tt)

labels = {'xdata': 'Time [s]', 'ydata': 'Signal amplitude [a.u.]',
          'xilt': r'$T_2$ [s]', 'yres': 'Residuals [a.u.]', 'titulo': 'CPMG'}

ilt = ILT(alpha=alpha, rango=(1e-4, 1e0), kernel=T2kernel, Nilt=Npts_L,
          figure=2, savepath=savepath, labels=labels, yscale='log')

ilt_distributions = []

for ii, expn in enumerate(expns):
    print(f"expn: {expn}")
    ydata = cpmgs[ii] / cpmgs[ii][0]
    ydata -= np.mean(ydata[-100:])
    ilt.DoTheStuff(ydata, times[ii], muestra=samples[ii])
    ilt_distributions.append(ilt.yilt)

ilt.legend()

ilt_distributions = np.asarray(ilt_distributions)
T2_axis = ilt.xilt

# Guardar ILT
if save:
    np.savez(f"{savepath}ILT_vs_experiment.npz",
             expns=expns, T2=T2_axis, distributions=ilt_distributions)

# ============================================================
# BIEXPONENTIAL FIT
# ============================================================

def biexp_offset(t, A1, T2a, A2, T2b, y0):
    return A1*np.exp(-t/T2a) + A2*np.exp(-t/T2b) + y0

T2short, T2long, Pshort, Plong, offsets = [], [], [], [], []

for ii in range(len(cpmgs)):
    t = times[ii] * 1000
    ydata = cpmgs[ii] / cpmgs[ii][0]
    ydata -= np.mean(ydata[-100:])
    p0 = [0.5, t[int(len(t)*0.3)], 0.5, t[int(len(t)*0.7)], 0]
    bounds = ([0, 0, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])

    try:
        popt, _ = curve_fit(biexp_offset, t, ydata, p0=p0, bounds=bounds, maxfev=50000)
        A1, T2a, A2, T2b, y0 = popt

        if T2a < T2b:
            ts, tl, As, Al = T2a, T2b, A1, A2
        else:
            ts, tl, As, Al = T2b, T2a, A2, A1

        total = As + Al
        T2short.append(ts)
        T2long.append(tl)
        Pshort.append(As / total)
        Plong.append(Al / total)
        offsets.append(y0)

    except RuntimeError:
        T2short.append(np.nan)
        T2long.append(np.nan)
        Pshort.append(np.nan)
        Plong.append(np.nan)
        offsets.append(np.nan)

T2short, T2long = np.asarray(T2short), np.asarray(T2long)
Pshort, Plong = np.asarray(Pshort), np.asarray(Plong)

# ============================================================
# T2 VS EXPERIMENT
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(expns, T2short, 'o-', label=r'$T_{2,short}$')
ax.plot(expns, T2long, 's-', label=r'$T_{2,long}$')
ax.set_xlabel("Experiment number")
ax.set_ylabel(r"$T_2$ [ms]")
ax.set_yscale('log')
ax.grid(alpha=0.3, which='both')
ax.legend()
ax.set_title("Biexponential T₂ components")
plt.tight_layout()
plt.show()

# ============================================================
# WEIGHTS VS EXPERIMENT
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(expns, Pshort*100, 'o-', label=r'$T_{2,short}$')
ax.plot(expns, Plong*100, 's-', label=r'$T_{2,long}$')
ax.set_xlabel("Experiment number")
ax.set_ylabel("Weight [%]")
ax.set_ylim(0, 100)
ax.grid(alpha=0.3)
ax.legend()
ax.set_title("Biexponential component weights")
plt.tight_layout()
plt.show()

#%%============================================================
# 2D MAP: ILT VS EXPERIMENT
# ============================================================

T2_ms = T2_axis * 1000

fig, ax = plt.subplots(figsize=(8, 6))
pcm = ax.pcolormesh(T2_ms, expns, ilt_distributions, shading='auto', vmax=0.05)
ax.set_xscale('log')
ax.set_xlabel(r'$T_2$ [ms]')
ax.set_ylabel('Experiment number')
ax.set_title(r'ILT $T_2$ map')
cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label('ILT amplitude')
plt.tight_layout()
plt.show()

# %%
