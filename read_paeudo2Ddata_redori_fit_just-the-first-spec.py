# -*- coding: utf-8 -*-
"""
Generalized N-peak Voigt fit with lmfit
Includes translucent components and residuals subplot
Saves amplitudes and errors
"""
import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import VoigtModel
from Datos import *

###############################################################################
# ----------------------- USER PARAMETERS ------------------------------------
###############################################################################

# expn = 103
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
# plot_individual_pairs = True
# plotRange = [-66, -76]

# slice_number = 0  # Which slice to fit

# # Define N peaks
# centers  = [-69.0, -71.0]      # ppm
# sigmas   = [0.3, 0.3]
# gammas   = [0.5, 0.5]
# amps_guess = [100, 50]         # initial guess for amplitudes

# # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## -TFSI in PEO-PTT
# expn = 103
# path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
# plot_individual_pairs = True
# plotRange = [-88, -94]
# region = "TFSI"
# # -------- Define peaks here (GENERALIZABLE) --------
# centers  = [-91.29, -92]      # ppm
# sigmas   = [0.0, 1]
# gammas   = [0.215, 0.5]
# amps_guess = [1200, 500]         # initial guess for amplitudes

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## -TFSI in PEO-PTT
expn = 91
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2026-02-07_PEO-PTT_solid-electrolyte/"
plot_individual_pairs = True
plotRange = [5, -7.5]
region = "TFSI"
# -------- Define peaks here (GENERALIZABLE) --------
centers  = [-1.61, -1.66]      # ppm
sigmas   = [0, 0.156]
gammas   = [0.086, 2]
amps_guess = [200, 200]




###############################################################################
# ----------------------- READ 2D DATA ---------------------------------------
###############################################################################
Npeaks = len(centers)
path_2D = f"{path}/{expn}/"
datos = DatosProcesados2D(path_2D, read_pp=False)
datos.espectro.ppmSelect(plotRange)

ppmAxis = datos.espectro.ppmAxis
spec = datos.espectro.real
slice_number = 1
ydata = spec[slice_number, :]
xdata = ppmAxis

###############################################################################
# ----------------------- BUILD N-PEAK MODEL ---------------------------------
###############################################################################

model = None
for i in range(Npeaks):
    prefix = f"p{i+1}_"
    peak = VoigtModel(prefix=prefix)
    model = peak if model is None else model + peak

###############################################################################
# ----------------------- SET PARAMETERS -------------------------------------
###############################################################################

params = model.make_params()

for i in range(Npeaks):
    prefix = f"p{i+1}_"
    params[f"{prefix}amplitude"].set(value=amps_guess[i], min=0)
    params[f"{prefix}center"].set(value=centers[i], vary=True)
    params[f"{prefix}sigma"].set(value=sigmas[i], vary=True)
    params[f"{prefix}gamma"].set(value=gammas[i], vary=True)

# TMP overwriting:
print("WARNING!!! FIXING SOME PEAKS!")
params[f"p1_center"].set(value=centers[0], vary=False)
params[f"p1_sigma"].set(value=sigmas[0], vary=False)
params[f"p1_gamma"].set(value=gammas[0], vary=False)


###############################################################################
# ----------------------- FIT ------------------------------------------------
###############################################################################

result = model.fit(ydata, params, x=xdata)

# Update guesses (optional)
for i in range(Npeaks):
    amps_guess[i] = result.params[f'p{i+1}_amplitude'].value

# Extract results
amps = []
amps_err = []
centers_fit = []
sigmas_fit = []
gammas_fit = []

for i in range(Npeaks):
    prefix = f"p{i+1}_"
    amp = result.params[f'{prefix}amplitude'].value
    err = result.params[f'{prefix}amplitude'].stderr or 0
    center = result.params[f'{prefix}center'].value
    sigma  = result.params[f'{prefix}sigma'].value
    gamma  = result.params[f'{prefix}gamma'].value

    amps.append(amp)
    amps_err.append(err*1.96)  # 95% confidence
    centers_fit.append(center)
    sigmas_fit.append(sigma)
    gammas_fit.append(gamma)

Stotal = np.sum(np.array(amps))
Stotal_err = np.sqrt(np.sum(np.array(amps_err)/1.96)**2)*1.96

###############################################################################
# ----------------------- PRINT RESULTS --------------------------------------
###############################################################################

print("=== Fit results ===")
for i in range(Npeaks):
    print(f"Peak {i+1}: {centers_fit[i]:.2f} ppm")
    print(f"  Amplitude = {amps[i]:.2f} ± {amps_err[i]:.2f}")
    print(f"  Sigma     = {sigmas_fit[i]:.3f}")
    print(f"  Gamma     = {gammas_fit[i]:.3f}")
print(f"Total signal: {Stotal:.2f} ± {Stotal_err:.2f}")

###############################################################################
# ----------------------- PLOT WITH COMPONENTS AND RESIDUALS -----------------
###############################################################################

fig, (ax1, ax2) = plt.subplots(
    2, 1, gridspec_kw={"height_ratios": [3, 1]}, figsize=(8, 5), sharex=True
)

# Data and fit
ax1.plot(xdata, ydata, 'k', label='Data')
ax1.plot(xdata, result.best_fit, 'r', label='Fit')

# Component peaks (translucent)
comps = result.eval_components(x=xdata)
for i in range(Npeaks):
    ax1.fill_between(
        xdata,
        comps[f'p{i+1}_'],
        alpha=0.3,
        label=f'Peak {i+1} ({centers_fit[i]:.2f} ppm)'
    )

ax1.set_ylabel("Intensity")
ax1.invert_xaxis()
ax1.legend(fontsize=8)
ax1.grid(True)

# Residuals
residuals = ydata - result.best_fit
ax2.plot(xdata, residuals, 'k')
ax2.axhline(0, color='gray', lw=1)
ax2.set_xlabel("ppm")
ax2.set_ylabel("Residuals")
ax2.grid(True)

plt.tight_layout()
plt.show()


