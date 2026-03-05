"""
Docstring for Exchange_Ramon_multi-fit
3-k exchange version




no anda!!!!!!!!!!!!!!!
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, report_fit
from Datos import *

#%% =========================
# Experimentos y temperaturas
# =========================
selected_experiments = [
    '13','22','32','42','52','62','72','82','92','102',
    '202','302','402','502','602','702','802','902'
]

temperatures = np.array([
    248.15,258.15,268.15,278.15,288.15,298.15,
    308.15,318.15,328.15,338.15,348.15,
    358.15,368.15,378.15,388.15,398.15,
    408.15,418.15
])

selected_experiments = selected_experiments[::10]
temperatures = temperatures[::10]

N = len(selected_experiments)

path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\4.0mm_6Li_Li3P_BM30_VT_NMR_10thJune2025/"
ppmRange = [9, 4]

# =========================
# Sitios en exchange
# =========================
nu01, nu02 = 533.9382, 333.12
fwhm1_init, fwhm2_init = 50, 50 #59.51, 59.97
fwhm_vary = True

p1 = 1 / (1 + 3.3566)
p2 = 1 - p1

# =========================
# Cargar espectros
# =========================
x_list, y_list, ppm_list = [], [], []

for expn in selected_experiments:
    datos = DatosProcesados(f'{path}{expn}/')
    datos.espectro.ppmSelect(ppmRange)

    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    ppm_list.append(ppmAxis)
    x_list.append(ppmAxis * datos.acqus.SFO1)
    y_list.append(spec)

#%% =========================
# Norrisian
# =========================
def norrisian(nu, k, amp, nu0s, fwhms, populations):

    Fn = []
    for n in range(len(nu0s)):
        Fn.append(
            populations[n] /
            (1j*(nu - nu0s[n]) - fwhms[n]/2 - k)
        )

    sumFn = np.sum(Fn, axis=0)
    Lex = -1j * sumFn
    Lex /= (1 + k * sumFn)
    Lex *= amp

    return np.imag(Lex)

#%% =========================
# Multi-k Norrisian (GENERAL)
# =========================
def multi_norrisian(nu, ks, Amp, weights, nu0s, fwhms, populations):

    Lex_tot = np.zeros_like(nu, dtype=float)

    for k, w in zip(ks, weights):
        Lex_tot += norrisian(
            nu,
            k,
            Amp * w,
            nu0s,
            fwhms,
            populations
        )

    return Lex_tot

#%% =========================
# Modelo por espectro
# =========================
def model_i(x, params, i):

    # Global fixed
    nu01 = params['nu01'].value
    nu02 = params['nu02'].value
    p1 = params['p1'].value
    p2 = params['p2'].value

    nu0s = [nu01, nu02]
    populations = [p1, p2]

    # global variable
    w1 = params[f'w1'].value
    w2 = params[f'w2'].value
    w3 = params[f'w3'].value
    weights = [w1, w2, w3]

    # Local
    fwhm1 = params[f'fwhm1_{i}'].value
    fwhm2 = params[f'fwhm2_{i}'].value
    fwhms = [fwhm1, fwhm2]

    Amp = params[f'Amp_{i}'].value

    k1 = params[f'k1_{i}'].value
    k2 = params[f'k2_{i}'].value
    k3 = params[f'k3_{i}'].value
    ks = [k1, k2, k3]

    # w1 = params[f'w1_{i}'].value
    # w2 = params[f'w2_{i}'].value
    # w3 = params[f'w3_{i}'].value
    # weights = [w1, w2, w3]

    # Impurities
    amp3 = params[f'amp3_{i}'].value
    amp4 = params[f'amp4_{i}'].value

    lor3 = norrisian(x, 0, amp3, [260.08], [45.25], [1])
    lor4 = norrisian(x, 0, amp4, [207.36], [61.76], [1])

    model = multi_norrisian(
        x,
        ks,
        Amp,
        weights,
        nu0s,
        fwhms,
        populations
    )

    return model + lor3 + lor4

#%% =========================
# Residuos globales
# =========================
def residuals(params, x_list, y_list):

    res = []
    for i, (x, y) in enumerate(zip(x_list, y_list)):
        res.append(y - model_i(x, params, i))

    return np.concatenate(res)

#%% =========================
# Parámetros lmfit
# =========================
params = Parameters()

# Global fixed
params.add('nu01', value=nu01, vary=False)
params.add('nu02', value=nu02, vary=False)
params.add('p1', value=p1, vary=False)
params.add('p2', value=p2, vary=False)


# Weights (sum = 1)
params.add(f'w1', value=0.4, min=0, max=1)
params.add(f'w2', value=0.3, min=0, max=1)
params.add(f'w3', expr=f'1 - w1 - w2', min=0, max=1)


# Local parameters
for i in range(N):

    params.add(f'fwhm1_{i}', value=fwhm1_init, min=fwhm1_init*0.9, max=fwhm1_init*1.1, vary=fwhm_vary)
    params.add(f'fwhm2_{i}', value=fwhm2_init, min=fwhm2_init*0.9, max=fwhm2_init*1.1, vary=fwhm_vary)

    params.add(f'Amp_{i}', value=1e3, min=0)

    # 3 exchange rates
    params.add(f'k1_{i}', value=0, min=0, max=10, vary=False)
    params.add(f'k2_{i}', value=100, min=0, max=500)
    params.add(f'k3_{i}', value=350, min=0, max=10000)

    # # Weights (sum = 1)
    # params.add(f'w1_{i}', value=0.4, min=0, max=1)
    # params.add(f'w2_{i}', value=0.3, min=0, max=1)
    # params.add(f'w3_{i}', expr=f'1 - w1_{i} - w2_{i}', min=0, max=1)

    # Impurities
    params.add(f'amp3_{i}', value=0, min=0, vary=False)
    params.add(f'amp4_{i}', value=0, min=0, vary=False)

#%% =========================
# Ajuste global
# =========================
result = minimize(
    residuals,
    params,
    args=(x_list, y_list),
    method='leastsq'
)

print("=== Ajuste global ===")
report_fit(result)

#%% =========================
# Graficar ajustes
# =========================
for i, (x, y) in enumerate(zip(x_list, y_list)):

    p = result.params

    k1 = p[f'k1_{i}'].value
    k2 = p[f'k2_{i}'].value
    k3 = p[f'k3_{i}'].value

    # w1 = p[f'w1_{i}'].value
    # w2 = p[f'w2_{i}'].value
    # w3 = p[f'w3_{i}'].value
    w1 = p[f'w1'].value
    w2 = p[f'w2'].value
    w3 = p[f'w3'].value

    Amp = p[f'Amp_{i}'].value
    fwhm1 = p[f'fwhm1_{i}'].value
    fwhm2 = p[f'fwhm2_{i}'].value

    nu0s = [p['nu01'].value, p['nu02'].value]
    pops = [p['p1'].value, p['p2'].value]
    fwhms = [fwhm1, fwhm2]

    comp1 = norrisian(x, k1, Amp*w1, nu0s, fwhms, pops)
    comp2 = norrisian(x, k2, Amp*w2, nu0s, fwhms, pops)
    comp3 = norrisian(x, k3, Amp*w3, nu0s, fwhms, pops)

    y_fit = comp1 + comp2 + comp3

    plt.figure(figsize=(8,5))
    plt.plot(x, y, '.', color='black', label='Data')
    plt.plot(x, y_fit, 'r-', label='Total fit')
    plt.plot(x, comp1, '--', label=f'k1 = {k1:.1f}')
    plt.plot(x, comp2, '--', label=f'k2 = {k2:.1f}')
    plt.plot(x, comp3, '--', label=f'k3 = {k3:.1f}')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Intensity')
    plt.title(f'T = {temperatures[i]-273.15:.1f} °C')
    plt.legend()
    plt.show()

#%% =========================
# k vs 1/T
# =========================
Tinv = 1 / temperatures

k1s = np.array([result.params[f'k1_{i}'].value for i in range(N)])
k2s = np.array([result.params[f'k2_{i}'].value for i in range(N)])
k3s = np.array([result.params[f'k3_{i}'].value for i in range(N)])

plt.figure(figsize=(8,5))
plt.plot(Tinv, k1s, 'o-', label='k1')
plt.plot(Tinv, k2s, 'o-', label='k2')
plt.plot(Tinv, k3s, 'o-', label='k3')
plt.yscale('log')
plt.xlabel('1 / T (1/K)')
plt.ylabel('k (s$^{-1}$)')
plt.legend()
plt.show()
