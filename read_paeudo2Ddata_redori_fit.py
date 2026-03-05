# -*- coding: utf-8 -*-
"""
REDOR analysis using lmfit
Generalized to N peaks
Stores amplitudes for later analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from Datos import *
from redor import *
from lmfit.models import VoigtModel
import pandas as pd
#################### functions ################################################

def label_curve(ax, x, y, label, idx, offset=(0, 0), **kwargs):
    """Annotate a curve at index idx (no rotation)."""
    x0 = x[idx] + offset[0]
    y0 = y[idx] + offset[1]

    ax.text(x0, y0, label,
            rotation=0,
            ha='center', va='center',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.5),
            **kwargs)
###############################################################################
# ----------------------- USER PARAMETERS ------------------------------------
###############################################################################

# # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## -CF3 in PEO-PTT
# expn = 103
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
# savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-02_500MHz_13C-CP\19F-7Li_REDOR/"
# save = True
# plot_individual_pairs = True
# plotRange = [-66, -76]
# max_recopl_time = 1.24 # ms
# spin_speed = 14000
# region = "_-CF3_of_PTT"
# nuclei = ['19F', '7Li']
# # -------- Define peaks here (GENERALIZABLE) --------
# centers = [-69.88, -70.50]
# sigmas  = [0.2965, 1.2004]
# gammas  = [0.4842, 0.3717]
# amp_guess = [40, 140]

#######- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## -TFSI in PEO-PTT
expn = 103
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-02_500MHz_13C-CP\19F-7Li_REDOR/"
save = True
plot_individual_pairs = True
plotRange = [-88, -94]
max_recopl_time = 10 # ms
spin_speed = 14000
nuclei = ['19F', '7Li']
region = "_TFSI"
# -------- Define peaks here (GENERALIZABLE) --------
centers = [-91.29, -91.32]
sigmas  = [0.0, 0.111]
gammas  = [0.189, 0.415]
amp_guess = [1282, 1647]


# # # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## 7Li{1H} redor
# expn = 91
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
# savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\PolymerElectrolyte\Analysis\2026-02_500MHz_13C-CP\7Li-1H_REDOR/"
# save = True
# plot_individual_pairs = True
# plotRange = [5, -7.5]
# max_recopl_time = 1.25 # ms
# spin_speed = 10000
# nuclei = ['7Li', '1H']
# region = ""
# # -------- Define peaks here (GENERALIZABLE) --------
# centers = [-1.61, -1.61]
# sigmas  = [0.0, 0.0]
# gammas  = [0.086, 1.729]
# amp_guess = [206, 406]

###############################################################################
# ----------------------- READ 2D DATA ---------------------------------------
###############################################################################
Npeaks = len(centers)
path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(path_2D, read_pp=False)
datos.espectro.ppmSelect(plotRange)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real

###############################################################################
# ----------------------- BUILD MODEL (N peaks) ------------------------------
###############################################################################

model = None
for i in range(Npeaks):
    prefix = f"p{i+1}_"
    peak = VoigtModel(prefix=prefix)
    model = peak if model is None else model + peak

###############################################################################
# ----------------------- FIT LOOP -------------------------------------------
###############################################################################

Signals_total = []
Signals_total_err = []

Signals_peaks = [[] for _ in range(Npeaks)]
Signals_peaks_err = [[] for _ in range(Npeaks)]

for kk in range(spec.shape[0]):

    ydata = spec[kk, :]
    xdata = ppmAxis

    params = model.make_params()

    # Assign parameters peak-by-peak
    for i in range(Npeaks):
        prefix = f"p{i+1}_"
        params[f'{prefix}amplitude'].set(value=amp_guess[i], min=0)
        params[f'{prefix}center'].set(value=centers[i], vary=False)
        params[f'{prefix}sigma'].set(value=sigmas[i], vary=False)
        params[f'{prefix}gamma'].set(value=gammas[i], vary=False)

    result = model.fit(ydata, params, x=xdata)

    # Update guesses iteratively
    for i in range(Npeaks):
        amp_guess[i] = result.params[f'p{i+1}_amplitude'].value

    # ----- Extract amplitudes -----
    amps = []
    errs = []

    for i in range(Npeaks):
        amp = result.params[f'p{i+1}_amplitude'].value
        err = result.params[f'p{i+1}_amplitude'].stderr or 0
        amps.append(amp)
        errs.append(err)

        Signals_peaks[i].append(amp)
        Signals_peaks_err[i].append(1.96 * err)

    Stotal = np.sum(amps)
    Stotal_err = 1.96 * np.sqrt(np.sum(np.array(errs)**2))

    Signals_total.append(Stotal)
    Signals_total_err.append(Stotal_err)

    # ---------- Plot with translucent components ----------
    if plot_individual_pairs:
        comps = result.eval_components(x=xdata)
        plt.figure()
        plt.plot(xdata, ydata, 'k', label='Data')
        plt.plot(xdata, result.best_fit, 'r', label='Fit')

        for i in range(Npeaks):
            plt.fill_between(
                xdata,
                comps[f'p{i+1}_'],
                alpha=0.3,
                label=f'Peak {i+1}'
            )

        plt.gca().invert_xaxis()
        plt.title(f"Slice {kk}")
        plt.legend()        

###############################################################################
# ----------------------- Convert to arrays ----------------------------------
###############################################################################

Signals_total = np.array(Signals_total)
Signals_total_err = np.array(Signals_total_err)

Signals_peaks = [np.array(p) for p in Signals_peaks]
Signals_peaks_err = [np.array(p) for p in Signals_peaks_err]

###############################################################################
# ----------------------- REDOR processing -----------------------------------
###############################################################################

S  = Signals_total[0::2]
S0 = Signals_total[1::2]
S_err  = Signals_total_err[0::2]
S0_err = Signals_total_err[1::2]

N = np.arange(1, len(S) + 1)
recopl_time_ms = (N / spin_speed) * 1000

###############################################################################
# ----------------------- SAVE DATA TABLE ------------------------------------
###############################################################################

# Separate S and S0 for total
S_total  = Signals_total[0::2]
S0_total = Signals_total[1::2]

S_total_err  = Signals_total_err[0::2]
S0_total_err = Signals_total_err[1::2]

data_dict = {
    "dephasing_time": recopl_time_ms,
    "Stotal": S_total,
    "Stotalerr": S_total_err,
    "S0total": S0_total,
    "S0totalerr": S0_total_err,
}

# -------- Individual peaks --------
for i in range(Npeaks):

    ppm_value = centers[i]

    S_i  = Signals_peaks[i][0::2]
    S0_i = Signals_peaks[i][1::2]

    S_i_err  = Signals_peaks_err[i][0::2]
    S0_i_err = Signals_peaks_err[i][1::2]

    data_dict[f"S{i+1}"] = S_i
    data_dict[f"S{i+1}err"] = S_i_err

    data_dict[f"S0{i+1}"] = S0_i
    data_dict[f"S0{i+1}err"] = S0_i_err

df = pd.DataFrame(data_dict)

if save:
    df.to_csv(f"{savepath}REDOR_amplitudes{region}.csv", index=False)

#%%##############################################################################
# ----------------------- REDOR TOTAL Y POR PICO --------------------------------
###############################################################################

# Filtrar hasta el tiempo máximo de recoupling
mask = recopl_time_ms <= max_recopl_time

# --- REDOR total ---
S_S0 = S[mask] / S0[mask]
S_S0_err = S_S0 * np.sqrt(
    (S_err[mask]/S[mask])**2 + (S0_err[mask]/S0[mask])**2
)

fig_redor, ax_redor = plt.subplots(figsize=(6,4))

#### Simulaciones teóricas (distancias)
for internuc_distance in np.arange(0.55, 0.8, 0.05):
    Nsim = np.arange(1, 1.2*N[-1])
    NTr, S_S0_sim = DeltaS_quadrupolar(internuc_distance, spin_speed, Nsim, nuclei=nuclei)
    xvals = NTr * 1000
    ax_redor.plot(xvals, S_S0_sim, 'k')
    label_curve(ax_redor, xvals, S_S0_sim,
                fr"{internuc_distance*10:.1f} $\AA$",
                idx=int(1.1*N[-1]), offset=(0,0.05), fontsize=8)

# Datos experimentales - total
ax_redor.errorbar(
    recopl_time_ms[mask], S_S0,
    yerr=S_S0_err,
    fmt='o', alpha=0.8, capsize=3,
    label='Total Experiment'
)

# # --- REDOR por pico ---
# for i in range(Npeaks):
#     S_i = Signals_peaks[i][0::2][mask]
#     S0_i = Signals_peaks[i][1::2][mask]
#     S_i_err = Signals_peaks_err[i][0::2][mask]
#     S0_i_err = Signals_peaks_err[i][1::2][mask]

#     S_i_S0 = S_i / S0_i
#     S_i_S0_err = S_i_S0 * np.sqrt(
#         (S_i_err / S_i)**2 + (S0_i_err / S0_i)**2
#     )

#     ax_redor.errorbar(
#         recopl_time_ms[mask], S_i_S0,
#         yerr=S_i_S0_err,
#         fmt='o', alpha=0.6, capsize=3,
#         label=f'Peak {i+1} ({centers[i]} ppm)'
#     )

ax_redor.set_xlabel(r"$\tau_{\text{REDOR}}$ (ms)")
ax_redor.set_ylabel("S/S0")
ax_redor.set_ylim(0, 1.1)
ax_redor.grid(True)
ax_redor.legend(fontsize=8)

#%%##############################################################################
# ----------------------- Plot S and S0 per peak --------------------------------
###############################################################################

fig, ax = plt.subplots(figsize=(7,5))

for i in range(Npeaks):
    # Filtrar hasta tiempo máximo
    mask = recopl_time_ms <= max_recopl_time

    S_i  = Signals_peaks[i][0::2][mask]
    S0_i = Signals_peaks[i][1::2][mask]
    S_i_err  = Signals_peaks_err[i][0::2][mask]
    S0_i_err = Signals_peaks_err[i][1::2][mask]

    # Plot S_i
    ax.errorbar(
        recopl_time_ms[mask],
        S_i,
        yerr=S_i_err,
        fmt='o-',
        alpha=0.7,
        capsize=3,
        label=f"S{i+1} ({centers[i]} ppm)"
    )

    # Plot S0_i
    ax.errorbar(
        recopl_time_ms[mask],
        S0_i,
        yerr=S0_i_err,
        fmt='s--',
        alpha=0.7,
        capsize=3,
        label=f"S0{i+1} ({centers[i]} ppm)"
    )

ax.set_xlabel(r"$\tau_{\text{REDOR}}$ (ms)")
ax.set_ylabel("Signal amplitude")
ax.set_title("S and S0 per peak vs Dephasing time")
ax.grid(True)
ax.legend(fontsize=8)
# ax.set_yscale('log')
# %%
