# -*- coding: utf-8 -*-
"""
Adaptado para usar PseudoVoigtFit en lugar de VoigtFit

@author: Santi


TRABAJAR EN LOS FITTINGS!!!!!
"""

import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rcParams['font.size'] = 12
import numpy as np
from Datos import *
import scipy.integrate as integrate
from VoigtFit import PseudoVoigtFit   # <--- cambio aquí
import pandas as pd


# directorio de datos
expns = np.arange(33,234)
expns_to_skip = np.arange(50, 366, 20)
expns = np.setdiff1d(expns, expns_to_skip)



absolute= False
autoph = False
path = rf"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\300old\2025-10-01_ccATMC_Rui-R3_NMC-Cu_PC/"
# directorio de guradado
savepath= r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\LiMetal\Rui\analysis\2025-10_R3/"
muestra = "7Li_cellR3-PCprotocol"

save = False
plotRange = [350, 150]
ppmRange = [300,200]

peaks = [245, 260]
range_of_peaks_to_save = [-50, -65]

# -------- parámetros iniciales 1 pico --------
m0_center = 240
m0_sigma = 10
m0_fraction = 0.5
m0_fwhm = 2.3548200*m0_sigma

# -------- parámetros iniciales 3 picos --------
m1_center, m2_center, m3_center = [240, 245, 250]
sigma = 6.5
m1_sigma,  m2_sigma,  m3_sigma  = sigma, sigma, sigma
m1_fraction, m2_fraction, m3_fraction = [0.5, 0.5, 0.5]
m1_fwhm, m2_fwhm, m3_fwhm = [2.3548200*m1_sigma, 2.3548200*m2_sigma, 2.3548200*m3_sigma]

# dataframe de resultados: ahora incluye fraction y el pico 3
vfit_results = pd.DataFrame(columns=[
    'expn', 'time',
    # 'm0_amplitude','m0_amplitude_stderr',
    # 'm0_center','m0_center_stderr',
    # 'm0_sigma','m0_sigma_stderr',
    # 'm0_fraction','m0_fraction_stderr',
    # 'm0_fwhm',
    'm1_amplitude','m1_amplitude_stderr',
    'm1_center','m1_center_stderr',
    'm1_sigma','m1_sigma_stderr',
    'm1_fraction','m1_fraction_stderr',
    'm1_fwhm',
    'm2_amplitude','m2_amplitude_stderr',
    'm2_center','m2_center_stderr',
    'm2_sigma','m2_sigma_stderr',
    'm2_fraction','m2_fraction_stderr',
    'm2_fwhm',
    'm3_amplitude','m3_amplitude_stderr',
    'm3_center','m3_center_stderr',
    'm3_sigma','m3_sigma_stderr',
    'm3_fraction','m3_fraction_stderr',
    'm3_fwhm',
    #'r-suqared_1peak', 
    'r-squared_3peaks'
])

colors = ['k', 'b', 'r', 'g', 'c', 'm', 'y']

fig_spec, ax_spec = plt.subplots(num=17856, nrows=1, figsize=(6, 4))
for jj, expn in enumerate(expns):
    print(f"Processing expn. {expn} ({jj+1}/{expns.size})")
    #=====================================================================
    # Ajuste de espectros 1D
    #=====================================================================
    datos = DatosProcesados(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real
    spec_time = datos.acqus.dic['DATE_START']

    spec1d_re = datos.espectro.real
    if jj==0:
        t_0 = spec_time
        spec_vs_t = np.zeros([spec.size, expns.size])
    spec_vs_t[:, jj] = spec1d_re
    ax_spec.plot(ppmAxis, spec1d_re)
    ax_spec.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]
    ax_spec.axvline(r1, color='k', linestyle='--')
    ax_spec.axvline(r2, color='k', linestyle='--')
    ax_spec.set_xlabel('chemical shift [ppm]')
    ax_spec.set_ylabel('Intensity [a.u.]')
    ax_spec.axhline(0, color='k')

    m0_amplitude = -np.trapezoid(spec1d_re, x=ppmAxis)
    if jj==0:
        m1_amplitude = m0_amplitude/3
        m2_amplitude = m0_amplitude/3
        m3_amplitude = m0_amplitude/3
    m0_height = spec1d_re.max()
    m1_height = m2_height = m3_height = m0_height/3


    # selecciono región de interés (ROI)
    semiancho = 60
    ppm_of_max = ppmAxis[spec1d_re == spec1d_re.max()][0]
    spec1d_ROI = spec1d_re[(ppmAxis > ppm_of_max-semiancho) & (ppmAxis < ppm_of_max+semiancho)]
    ppmAxis_ROI = ppmAxis[(ppmAxis > ppm_of_max-semiancho) & (ppmAxis < ppm_of_max+semiancho)]
    # -------------------------------------------
    # -------- ajuste con 3 PseudoVoigts -------- 
    vfit = PseudoVoigtFit(
        ppmAxis_ROI, spec1d_ROI,
        Npicos=3,
        ajustar=True,
        amplitude=[m1_amplitude, m2_amplitude, m3_amplitude],
        center=[m1_center, m2_center, m3_center],
        sigma=[m1_sigma, m2_sigma, m3_sigma],
        fraction=[m1_fraction, m2_fraction, m3_fraction],
        height=[m1_height, m2_height, m3_height],
        fwhm=[m1_fwhm, m2_fwhm, m3_fwhm],
        fijar = ["sigma"],
        bounds={"fraction": [0.2, 1],
                "m1_center": [m1_center-2.5, m1_center+2.5],
                "m2_center": [m2_center-2.5, m2_center+2.5],
                "m3_center": [m3_center-2.5, m3_center+2.5]}  # <--- aquí pongo límites en los centros
    )
    fig = vfit.plot_ajuste()
    fig.gca().set_title(f"expn. {expn} - 3 peaks PseudoVoigt fit")
    fig.savefig(f"{savepath}/3peaksFit_{muestra}_expn{expn}.png", dpi=300)
    plt.close(fig)
    ####redefino parámetros para la próxima vuelta
    # for k in [1,2,3]:
    #     globals()[f"m{k}_amplitude"] = vfit.params[f"m{k}_amplitude"].value
    #     globals()[f"m{k}_center"]    = vfit.params[f"m{k}_center"].value
    #     globals()[f"m{k}_sigma"]     = vfit.params[f"m{k}_sigma"].value
    #     globals()[f"m{k}_fraction"]  = vfit.params[f"m{k}_fraction"].value
    #     globals()[f"m{k}_fwhm"]      = 2.3548200 * globals()[f"m{k}_sigma"]
    #     globals()[f"m{k}_height"]    = vfit.params[f"m{k}_height"].value
    # -------------------------------------------

    # # --------------------------------------------
    # # -------- ajuste con 1 PseudoVoigt --------
    # # centro como promedio ponderado de los 3 picos
    # total_amp = sum([vfit.params[f"m{k}_amplitude"].value for k in [1,2,3]])
    # m0_center = sum([vfit.params[f"m{k}_amplitude"].value * vfit.params[f"m{k}_center"].value for k in [1,2,3]]) / total_amp

    # vfit0 = PseudoVoigtFit(
    #     ppmAxis_ROI, spec1d_ROI,
    #     Npicos=1,
    #     ajustar=False,
    #     amplitude=[m0_amplitude],
    #     center=[m0_center],
    #     sigma=[m0_sigma],
    #     fraction=[m0_fraction],
    #     height=[m0_height],
    #     fwhm=[m0_fwhm]
    # )
    # # fig0 = vfit0.plot_ajuste()   # <--- gris
    # # fig0.gca().set_title(f"expn. {expn} - 1 peak PseudoVoigt fit")
    # # --------------------------------------------

    # guardo resultados en dataframe
    df = pd.DataFrame({
        'expn': [expn],
        'time': spec_time,
        # 'm0_amplitude': [vfit0.params["m1_amplitude"].value], 
        # 'm0_amplitude_stderr': [vfit0.params["m1_amplitude"].stderr],
        # 'm0_center': [vfit0.params["m1_center"].value],
        # 'm0_center_stderr': [vfit0.params["m1_center"].stderr],
        # 'm0_sigma': [vfit0.params["m1_sigma"].value],
        # 'm0_sigma_stderr': [vfit0.params["m1_sigma"].stderr],
        # 'm0_fraction': [vfit0.params["m1_fraction"].value],
        # 'm0_fraction_stderr': [vfit0.params["m1_fraction"].stderr],
        # 'm0_fwhm': [vfit0.params["m1_fwhm"].value],
        'm1_amplitude': [vfit.params["m1_amplitude"].value],
        'm1_amplitude_stderr': [vfit.params["m1_amplitude"].stderr],
        'm1_center': [vfit.params["m1_center"].value],
        'm1_center_stderr': [vfit.params["m1_center"].stderr],
        'm1_sigma': [vfit.params["m1_sigma"].value],
        'm1_sigma_stderr': [vfit.params["m1_sigma"].stderr],
        'm1_fraction': [vfit.params["m1_fraction"].value],
        'm1_fraction_stderr': [vfit.params["m1_fraction"].stderr],
        'm1_fwhm': [vfit.params["m1_fwhm"].value],
        'm2_amplitude': [vfit.params["m2_amplitude"].value],
        'm2_amplitude_stderr': [vfit.params["m2_amplitude"].stderr],
        'm2_center': [vfit.params["m2_center"].value],
        'm2_center_stderr': [vfit.params["m2_center"].stderr],
        'm2_sigma': [vfit.params["m2_sigma"].value],
        'm2_sigma_stderr': [vfit.params["m2_sigma"].stderr],
        'm2_fraction': [vfit.params["m2_fraction"].value],
        'm2_fraction_stderr': [vfit.params["m2_fraction"].stderr],
        'm2_fwhm': [vfit.params["m2_fwhm"].value],
        'm3_amplitude': [vfit.params["m3_amplitude"].value],
        'm3_amplitude_stderr': [vfit.params["m3_amplitude"].stderr],
        'm3_center': [vfit.params["m3_center"].value],
        'm3_center_stderr': [vfit.params["m3_center"].stderr],
        'm3_sigma': [vfit.params["m3_sigma"].value],
        'm3_sigma_stderr': [vfit.params["m3_sigma"].stderr],
        'm3_fraction': [vfit.params["m3_fraction"].value],
        'm3_fraction_stderr': [vfit.params["m3_fraction"].stderr],
        'm3_fwhm': [vfit.params["m3_fwhm"].value],
        # 'r-suqared_1peak': [vfit0.ajuste.rsquared],
        'r-squared_3peaks': [vfit.ajuste.rsquared]
    })
    vfit_results = pd.concat([vfit_results, df])

vfit_results.to_csv(f"{savepath}/vfit_results_3peaks_{muestra}.csv", index=False)
#%% -------------
fig, ax = plt.subplots(num=1785731, figsize=(8, 3))
fig_int, ax_int = plt.subplots(num=382910, figsize=(8, 3))

time = vfit_results['time'] - vfit_results['time'].min()
time = time/3600  # hours

signals = vfit_results['m1_amplitude'] + vfit_results['m2_amplitude'] + vfit_results['m3_amplitude']

m1_center = vfit_results['m1_center']
m2_center = vfit_results['m2_center']
m3_center = vfit_results['m3_center']
m1_amplitude = vfit_results['m1_amplitude']
m2_amplitude = vfit_results['m2_amplitude']
m3_amplitude = vfit_results['m3_amplitude']

# ax_int.plot(time, vfit_results['m0_amplitude']/signals.max(), 'o-', color='gray', label='total signal: 1 peak fit')
ax_int.plot(time, m1_amplitude/signals.max(), 'ro', label='peak 1: 240 ppm')
ax_int.plot(time, m2_amplitude/signals.max(), 'bo', label='peak 2: 247 ppm')
ax_int.plot(time, m3_amplitude/signals.max(), 'go', label='peak 3: 256 ppm')
ax_int.plot(time, signals/signals.max(), 'ko', label='total signal')
ax_int.yaxis.set_major_locator(MultipleLocator(0.1))
ax_int.grid(axis='y')
ax_int.set_xlabel('Time [h]')
ax_int.set_ylabel('Normalized area')
ax_int.set_ylim(-0.05, 1.1)
ax_int.legend()

# ax.plot(time, vfit_results['m0_center'], 'o', color='gray', label='peak 0 (1 peak fit)')
ax.plot(time, m1_center, 'ro', label='peak 1: 237 ppm')
ax.plot(time, m2_center, 'bo', label='peak 2: 245 ppm')
ax.plot(time, m3_center, 'go', label='peak 3: 254 ppm')
ax.set_xlabel('Time [h]')
ax.set_ylabel(r'$\delta$ [ppm]')
ax.set_ylim(232, 258)
ax.legend()

#%%
ppm_mesh, tau_mesh = np.meshgrid(ppmAxis, time)
fig, ax = plt.subplots(num=17856522, figsize=(3, 6))
pcm = ax.pcolormesh(ppm_mesh, tau_mesh, spec_vs_t.T, shading='auto')
ax.set_ylabel('Time [h]')
ax.set_xlabel(r'$^7$Li Chemical shift [ppm]')
fig.colorbar(pcm, ax=ax, label='Intensity')
ax.set_xlim([260, 220])

#%%

fig_par, ax_par = plt.subplots(num=1785652313, figsize=(8, 3))
# ax_par.plot(time, vfit_results['r-suqared_1peak'], 'o-', color='gray', label='r-squared (1 peak fit)')
ax_par.plot(time, vfit_results['r-squared_3peaks'], 'ko-', label='r-squared (3 peaks fit)')
ax_par.set_ylim([0.98, 1.001])
ax_par.set_ylabel(r'$R^2$')
ax_par.set_xlabel('Time [h]')
ax_par.legend()

# %%

# fig_fwhm, ax_fwhm = plt.subplots(num=1785652313, figsize=(8, 3))
# fwhm_G = 2*np.sqrt(2*np.log(2))*vfit_results["m0_sigma"]
# fwhm_L = 2*vfit_results["m0_sigma"] 
# fwhm = fwhm_L/2 + np.sqrt(fwhm_L**2/4 + fwhm_G**2)  
# ax_fwhm.plot(time, fwhm, 'o-', color='gray', label='FWHM peak 0 (1 peak fit)')
# ax_fwhm.set_xlabel('Time [h]')
# ax_fwhm.set_ylabel('FWHM (ppm)')
# ax_fwhm.set_ylim([25, 50])
# ax_fwhm.legend()

# %%
