# -*- coding: utf-8 -*-
"""
Simultaneous pseudo-Voigt fit of multiple NMR spectra
Replaces individual fits with a global simultaneous WEIGHTED fit.
The weight for each spectrum is based on its integral (to balance spectra with different SNR), normalized by the mean integral.
Centers, FWHM and fractions are shared; amplitudes and offsets are per spectrum.
@author: Santi
"""

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Datos import *
import lmfit
from lmfit import Parameters
from lmfit.models import PseudoVoigtModel

#=====================================================================
# Directorios y parámetros
#=====================================================================


expns = [67,71,62,68,63,70,60,64,66,69]
expns = [70]
path  = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\400dnp\2025-11-03_3.2mm_Debashis-dendrites/"
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\DNP\Debashis\Analysis\2025-11_R6/"
muestra = ""
save = False

ppmRange  = [15, -15]   # fitting region
show_individual_fits = True
#=====================================================================
# Load spectra
#=====================================================================

x_list = []
y_list = []
integral_list = []
contact_times = []

for expn in expns:
    datos = DatosProcesados(f'{path}/{expn}/', read_pp=False)
    datos.espectro.ppmSelect(ppmRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    contact_time = datos.acqus.dic["P"][15]

    xdata = ppmAxis
    ydata = spec

    x_list.append(xdata)
    y_list.append(ydata)
    integral_list.append(np.abs(integrate.simpson(ydata, x=xdata)))
    contact_times.append(contact_time)

#=====================================================================
# Define residual function for global fitting
#=====================================================================

def Residuals_global(params, x_list, y_list, integral_list=None, n_peaks=2):
    """
    Residual function for simultaneous fitting of multiple spectra.
    Shared: centers, fwhm, fractions
    Per-spectrum: amplitudes and offsets
    """
    res = []
    n_spec = len(x_list)
    if integral_list is not None:
        ## Normalize integrals to create weights list
        integral_array = np.array(integral_list)
        weights = integral_array / np.mean(integral_array)
    else:
        weights = np.ones(n_spec)

    for s in range(n_spec):
        x = x_list[s]
        y = y_list[s]
        model_eval = np.zeros_like(x)

        for i in range(n_peaks):
            # per-spectrum amplitude
            amp = params[f'ampl_{s}_{i}'].value

            # shared parameters
            center = params[f'center_{i}'].value
            fwhm = params[f'fwhm_{i}'].value
            fraction = params[f'fraction_{i}'].value
            sigma = fwhm / 2
            gamma = fwhm / 2

            # construct lmfit Parameters for this component
            comp_params = Parameters()
            comp_params.add('amplitude', value=amp)
            comp_params.add('center', value=center)
            comp_params.add('sigma', value=sigma)
            comp_params.add('gamma', value=gamma)
            comp_params.add('fraction', value=fraction)

            model = PseudoVoigtModel()
            model_eval += model.eval(params=comp_params, x=x)

        # add per-spectrum offset
        offset = params[f'offset_{s}'].value
        
        res.append(np.sqrt(weights[s]) * (y - (model_eval + offset)))

    return np.concatenate(res)

#=====================================================================
# Initialize parameters
#=====================================================================

n_spec = len(x_list)
n_peaks = 2

# Initial guesses from previous A code
amplitude_guess = [ [y_list[s].max()/2, y_list[s].max()/2] for s in range(n_spec)]
center_guess = [4, -1]
fwhm_guess   = [5, 5]
fraction_guess = [0.5, 0.5]

params = Parameters()

# shared parameters
for i in range(n_peaks):
    params.add(f'center_{i}', value=center_guess[i], vary=True)
    params.add(f'fwhm_{i}', value=fwhm_guess[i], vary=True, min=0)
    params.add(f'fraction_{i}', value=fraction_guess[i], vary=True, min=0, max=1)

# per-spectrum amplitudes and offsets
for s in range(n_spec):
    for i in range(n_peaks):
        params.add(f'ampl_{s}_{i}', value=amplitude_guess[s][i], vary=True, min=0)
    params.add(f'offset_{s}', value=0.0, vary=True, min=0)

#=====================================================================
# Perform global minimization
#=====================================================================

minim = lmfit.Minimizer(Residuals_global, params, fcn_args=(x_list, y_list))
print("minimizando...")
result = minim.minimize()
lmfit.report_fit(result)

#=====================================================================
# Extract results per spectrum
#=====================================================================

results = []

for s in range(n_spec):
    x = x_list[s]
    y = y_list[s]
    yfit = np.zeros_like(x)
    ycomp_list = []

    for i in range(n_peaks):
        amp = result.params[f'ampl_{s}_{i}'].value
        center = result.params[f'center_{i}'].value
        fwhm = result.params[f'fwhm_{i}'].value
        fraction = result.params[f'fraction_{i}'].value
        sigma = fwhm / 2
        gamma = fwhm / 2

        # construct lmfit Parameters for this component
        comp_params = Parameters()
        comp_params.add('amplitude', value=amp)
        comp_params.add('center', value=center)
        comp_params.add('sigma', value=sigma)
        comp_params.add('gamma', value=gamma)
        comp_params.add('fraction', value=fraction)

        model = PseudoVoigtModel()
        ycomp = model.eval(params=comp_params, x=x)
        ycomp_list.append(ycomp)
        yfit += ycomp

    # add per-spectrum offset
    offset = result.params[f'offset_{s}'].value
    yfit += offset

    results.append({
        'expn': expns[s],
        'contact_time': contact_times[s],
        'm1_area': result.params[f'ampl_{s}_0'].value,
        'm2_area': result.params[f'ampl_{s}_1'].value,
        'm1_center': result.params[f'center_0'].value,
        'm2_center': result.params[f'center_1'].value,
        'm1_fwhm': result.params[f'fwhm_0'].value,
        'm2_fwhm': result.params[f'fwhm_1'].value,
        'm1_fraction': result.params[f'fraction_0'].value,
        'm2_fraction': result.params[f'fraction_1'].value,
        'offset': offset,
        'ydata': y,
        'ycomp1': ycomp_list[0],
        'ycomp2': ycomp_list[1],
        'yfit': yfit
    })

#=====================================================================
# Create DataFrame and sort
#=====================================================================

df_results = pd.DataFrame(results)
df_results.sort_values(by='contact_time', inplace=True, ignore_index=True)
print(df_results)

#%%=====================================================================
# Plot Areas vs Contact Time
#=====================================================================

colors = ['r', 'b']
fig_area, ax_area = plt.subplots()
ax_area.plot(df_results['contact_time'], df_results['m1_area'], 'o-', color=colors[0], label='Peak 1')
ax_area.plot(df_results['contact_time'], df_results['m2_area'], 'o-', color=colors[1], label='Peak 2')
ax_area.set_xlabel(r"Contact time ($\mu$s)")
ax_area.set_ylabel("Amplitude (a.u.)")
ax_area.legend()
ax_area.grid(True)

#%%=====================================================================
# Stack Plot of spectra with components
#=====================================================================
textshift = 30 # original: 30

colors = ['r', 'b']
fig_stack, ax_stack = plt.subplots(figsize=(3,8)) # original figsize=(3,8)
contact_times_array = df_results['contact_time'].values
ydata_array = np.array(df_results['ydata'].to_list())
ycomp1_array = np.array(df_results['ycomp1'].to_list())
ycomp2_array = np.array(df_results['ycomp2'].to_list())
yfit_array = np.array(df_results['yfit'].to_list())

# sort by contact_time
sort_idx = np.argsort(contact_times_array)
contact_times_sorted = contact_times_array[sort_idx]
ydata_sorted = ydata_array[sort_idx,:]
ycomp1_sorted = ycomp1_array[sort_idx,:]
ycomp2_sorted = ycomp2_array[sort_idx,:]
yfit_sorted = yfit_array[sort_idx,:]

for i, t in enumerate(contact_times_sorted):
    offset = 1
    n = np.max(ydata_sorted[i])
    ax_stack.plot(x_list[0], ydata_sorted[i]/n + i*offset, color='grey', lw=4)
    ax_stack.plot(x_list[0], yfit_sorted[i]/n + i*offset, color='k', lw=1)
    ax_stack.plot(x_list[0], ycomp1_sorted[i]/n + i*offset, color=colors[0], lw=0.8)
    ax_stack.plot(x_list[0], ycomp2_sorted[i]/n + i*offset, color=colors[1], lw=0.8)
    text_x = np.max(x_list[0]) - textshift
    ax_stack.text(text_x, i*offset + 0.7, f"{t:.0f} µs", va='center', ha='right', fontsize=12)

ax_stack.set_xlabel("Chemical Shift (ppm)")
ax_stack.set_ylabel("Stacked Spectra")
ax_stack.set_xlim(np.max(x_list[0]), np.min(x_list[0]))
ax_stack.grid(True)


#%%=====================================================================
# Gráfico de espectros superpuestos normalizados (solo datos)
# con colormap viridis y colorbar discreta con labels por contact time
#=====================================================================
import matplotlib.colors as mcolors

fig_super, ax_super = plt.subplots(figsize=(6,4), num=43904812)

# Extraer datos
contact_times_array = df_results["contact_time"].values
ydata_array = np.array(df_results["ydata"].to_list())

# Normalizar cada espectro por su máximo
ydata_norm = ydata_array / ydata_array.max(axis=1)[:, None]

# Generar una lista de colores, uno por espectro
num_spectra = len(contact_times_array)
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 1, num_spectra))

# Graficar espectros
offset = 0.2
for i, (t, color) in enumerate(zip(contact_times_array, colors)):
    ax_super.plot(x_list[0], ydata_norm[i]+offset*i, color=color, lw=3, alpha=0.9)
# Formato del gráfico
ax_super.set_xlim(np.max(x_list[0]), np.min(x_list[0]))
ax_super.set_xlabel("Chemical Shift (ppm)")
ax_super.set_ylabel("Normalized Intensity (a.u.)")
ax_super.axvline(df_results['m1_center'][0], color='k', lw=1, ls='--')
ax_super.axvline(df_results['m2_center'][0], color='k', lw=1, ls='--')
ax_super.grid(True)

# Crear colorbar discreta con un color por espectro
cmap_discrete = mcolors.ListedColormap(colors)
bounds = np.arange(num_spectra + 1)
norm = mcolors.BoundaryNorm(bounds, cmap_discrete.N)
sm = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm)
sm.set_array([])

# Agregar colorbar con etiquetas correspondientes a contact times
cbar = fig_super.colorbar(sm, ax=ax_super, ticks=bounds[:-1] + 0.5)
cbar.ax.set_yticklabels([f"{int(t)}" for t in contact_times_array])
cbar.set_label("Contact time (µs)")


#=====================================================================
# Save figures and results
#=====================================================================

if save:
    fig_stack.savefig(f"{savepath}/{muestra}_StackSpectra.png")
    fig_area.savefig(f"{savepath}/{muestra}_Areas_vs_ContactTime.png")
    df_results.to_csv(f"{savepath}/{muestra}_VoigtFit_Global_results.csv", index=False)
    cp_data = np.array([df_results['contact_time'], df_results['m1_area'], df_results['m2_area']]).T
    np.savetxt(f"{savepath}/{muestra}_VoigtFit_Global_cp_data.txt", cp_data)

# %%
