"""
Docstring for Exchange_Ramon_multi-fit
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize, report_fit
from Datos import *



# =========================
# Experimentos y temperaturas
# =========================
selected_experiments = ['13','22','32','42','52','62','72','82','92','102',
                        '202','302','402','502','602','702','802','902']


temperatures = np.array([248.15,258.15,268.15,278.15,288.15,298.15,
                          308.15,318.15,328.15,338.15,348.15,
                          358.15,368.15,378.15,388.15,398.15,
                          408.15,418.15])
savepath = r"C:\Users\Santi\OneDrive - University of Cambridge\Projects\DNP\Ramon\Exchange\spec/"

selected_experiments = selected_experiments[::10]
temperatures = temperatures[::10]

N = len(selected_experiments)
path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\4.0mm_6Li_Li3P_BM30_VT_NMR_10thJune2025/"
ppmRange = [9,4]

# =========================
# Sitios en exchange
# =========================
nu01, nu02  = 533.9382, 333.12
fwhm1, fwhm2 = 59.51, 59.97
fwhm_vary = True

p1 = 1 / (1 + 2)
p1 = 1 / (1 + 3.3566)
p2 = 1 - p1


# Inicializo
x_list = []
y_list = []
ppm_list = []
for expn in selected_experiments:
    datos = DatosProcesados(f'{path}{expn}/')
    datos.espectro.ppmSelect(ppmRange)
    ppmAxis = datos.espectro.ppmAxis
    spec = datos.espectro.real

    ppm_list.append(ppmAxis)    
    x_list.append(ppmAxis*datos.acqus.SFO1)
    y_list.append(spec)
    np.savetxt(f'{savepath}spec_T-{temperatures[selected_experiments.index(expn)]:.2f}K.dat', np.column_stack((ppmAxis, spec, x_list[-1])), header='ppm\tintensity\tFrequency(Hz)')



def norrisian(nu, k, amp, nu0s, fwhms, populations):
    """
    nu : frequency axis, Hz
    nu0s, fwhms, populations : lists
    """
    Fn = []
    for n in range(len(nu0s)):
        nu0 = nu0s[n]
        fwhm =fwhms[n]
        Pn = populations[n]
        Fn.append( Pn / (1j*(nu-nu0) - fwhm/2 - k) )
    sumFn = np.sum(Fn, axis=0)
    Lex = -1j * sumFn
    Lex /= (1 + k * sumFn)
    Lex *= amp
    return np.imag(Lex)

def bi_norrisian(nu, ka, kb, Amp, ratio, fwhm1, fwhm2, nu01, nu02, p1, p2):
    nu0s = [nu01, nu02]
    populations = [p1,p2]
    fwhms = [fwhm1, fwhm2]

    Lex_a = norrisian(nu, ka, Amp*ratio, nu0s, fwhms, populations)
    Lex_b = norrisian(nu, kb, Amp*(1-ratio), nu0s, fwhms, populations)


    return Lex_a + Lex_b




# ----------------------------
# 2. Definir modelo para cada espectro
# ----------------------------

def model_i(x, params, i):
    """
    Modelo para el espectro i.
    Usa parámetros globales y locales:
      - center: global ajustable
      - sigma: global fijo
      - amp_i: local ajustable
      - offset_i: local ajustable
    """
    # global fixed parameters:
    p1 = params['p1'].value
    p2 = params['p2'].value
    nu01 = params['nu01'].value
    nu02 = params['nu02'].value
    # global variable parameters:    
    ratio = params['ratio'].value
    # local parameters:
    fwhm1 = params[f'fwhm1_{i}'].value
    fwhm2 = params[f'fwhm2_{i}'].value
    Amp  = params[f'Amp_{i}'].value
    ka =  params[f'ka_{i}'].value
    kb =  params[f'kb_{i}'].value

    ## Impurities:
    amp3 = params[f'amp3_{i}'].value
    amp4 = params[f'amp4_{i}'].value
    
    lorentzian3 = norrisian(x, 0, amp3, [260.08], [45.25], [1])
    lorentzian4 = norrisian(x, 0, amp4, [207.36], [61.76], [1])

    model_i = bi_norrisian(x, ka, kb, Amp, ratio, fwhm1, fwhm2, nu01, nu02, p1, p2)
    model_i += lorentzian3
    model_i += lorentzian4   
    
    return model_i

# ----------------------------
# 3. Función de residuos global
# ----------------------------

def residuals(params, x_list, y_list):
    """
    Devuelve los residuos concatenados de todos los espectros.
    lmfit minimizará la suma de estos residuos al cuadrado.
    """
    res = []
    for i, (x, y) in enumerate(zip(x_list, y_list)):
        y_model = model_i(x, params, i)
        res.append(y - y_model)
    return np.concatenate(res)

# ----------------------------
# 4. Definir parámetros
# ----------------------------

params = Parameters()
# Parámetros globales fijos
params.add('p1', value=p1, vary=False)  # fijo durante el ajuste
params.add('p2', value=p2, vary=False)  # fijo durante el ajuste
params.add('nu01', value=nu01, vary=False)  # fijo durante el ajuste
params.add('nu02', value=nu02, vary=False)  # fijo durante el ajuste

# Parámetros globales ajustables
params.add('ratio', value=0.3, min=0.0, max=1.0, vary=False)  # se ajustará para todos los espectros

# Parámetros locales ajustables
for i in range(N):
    params.add(f'fwhm1_{i}', value=fwhm1, min=10, max=80, vary=fwhm_vary)
    params.add(f'fwhm2_{i}', value=fwhm2, min=10, max=80, vary=fwhm_vary)         
    params.add(f'Amp_{i}', value=1e3, min=0)         
    params.add(f'ka_{i}', value=300, min=0)         
    params.add(f'kb_{i}', value=1, min=0)
    # params.add(f'amp3_{i}', value=1600, min=0)
    # params.add(f'amp4_{i}', value=1400, min=0)
    params.add(f'amp3_{i}', value=0, min=0, vary=False)
    params.add(f'amp4_{i}', value=0, min=0, vary=False)


# ----------------------------
# 5. Ajuste global con lmfit
# ----------------------------

result = minimize(residuals, params, args=(x_list, y_list), method='leastsq')

# ----------------------------
# 6. Resultados
# ----------------------------

print("=== Ajuste global de N espectros ===")
report_fit(result)  # imprime todos los parámetros ajustados



##############################################
##################################################
##################################################
##################################################
###################################################

#%%----------------------------
# Graficar componentes del ajuste para cada espectro
# ----------------------------

colors = plt.cm.viridis(np.linspace(0,1,N))

for i, (x, y) in enumerate(zip(x_list, y_list)):
    # Recuperar parámetros optimizados
    params_i = result.params
    ka = params_i[f'ka_{i}'].value
    kb = params_i[f'kb_{i}'].value
    Amp = params_i[f'Amp_{i}'].value
    ratio = params_i['ratio'].value
    fwhm1 = params_i[f'fwhm1_{i}'].value
    fwhm2 = params_i[f'fwhm2_{i}'].value
    nu01 = params_i['nu01'].value
    nu02 = params_i['nu02'].value
    p1 = params_i['p1'].value
    p2 = params_i['p2'].value
    amp3 = params_i[f'amp3_{i}'].value
    amp4 = params_i[f'amp4_{i}'].value

    # Componentes de la bi-norrisian
    bi_a = norrisian(x, ka, Amp*ratio, [nu01, nu02], [fwhm1, fwhm2], [p1,p2])
    bi_b = norrisian(x, kb, Amp*(1-ratio), [nu01, nu02], [fwhm1, fwhm2], [p1,p2])

    # Impurezas
    lor3 = norrisian(x, 0, amp3, [260.08], [45.25], [1])
    lor4 = norrisian(x, 0, amp4, [207.36], [61.76], [1])

    # Modelo total
    y_fit = bi_a + bi_b + lor3 + lor4

    # Graficar
    plt.figure(figsize=(8,5))
    plt.plot(x, y, '.', color='black', label='Data')
    plt.plot(x, y_fit, '-', color='red', label='Total fit')
    plt.plot(x, bi_a, '--', color='blue', label=f'ka component (ka={ka:.1f})')
    plt.plot(x, bi_b, '--', color='green', label=f'kb component (kb={kb:.1f})')
    plt.plot(x, lor3, ':', color='purple', label='Lorentzian3')
    plt.plot(x, lor4, ':', color='orange', label='Lorentzian4')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Intensity')
    plt.title(f'Spectrum {i} - Temperature: {temperatures[i]-273.15:.1f} °C')
    plt.legend()
    plt.show()

# %%
fig, ax = plt.subplots(figsize=(8,5))
kas = np.array([result.params[f'ka_{ii}'].value for ii in range(N)])
kbs = np.array([result.params[f'kb_{ii}'].value for ii in range(N)])
Tinv = 1 / temperatures
ax.plot(Tinv, kas, 'o-', label="ka")
ax.plot(Tinv, kbs, 'o-', label="kb")
ax.set_xlabel('1 / T (1/K)')
ax.set_ylabel('k (s$^{-1}$)')
ax.legend()
ax.set_yscale('log')
# ax.set_ylim(1e0, 1e4)


# %%
fig, ax = plt.subplots(figsize=(8,5))
fwhm1 = np.array([result.params[f'fwhm1_{ii}'].value for ii in range(N)])
fwhm2 = np.array([result.params[f'fwhm2_{ii}'].value for ii in range(N)])
Tinv = 1 / temperatures
ax.plot(Tinv, 2/fwhm1, 'o-', label="T2_1 (s)")
ax.plot(Tinv, 2/fwhm2, 'o-', label="T2_2 (s)")
ax.set_xlabel('1 / T (1/K)')
ax.set_ylabel('T_2 (s)')
ax.legend()
ax.set_yscale('log')

# %%
