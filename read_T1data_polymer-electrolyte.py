import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from scipy.optimize import curve_fit

#========================= INPUT =======================================
fit_biexponential = False
####### 1H
# nexp_base = 31
# plotRange = [100, -100]
####### 19F
# nexp_base = 33
# plotRange = [-20, -150]
####### 7Li
nexp_base = 35
plotRange = [100, -100]

T_nominal = np.arange(-20, 41, 10)
T_real = [-17.8, -9.9, 0, 10, 20, 30, 40]

path = r"C:\\Users\\Santi\\OneDrive - University of Cambridge\\NMRdata\\500\\2025-07-16_PEO-solid-electrolyte_VT/"
ppmIntegrationWidth = 30

#========================= INIT ========================================
expns = [nexp_base + 10 * i for i in range(len(T_nominal))]
T1mono_list, T1mono_err_list = [], []
T1short_list, T1long_list, P1short_list = [], [], []

#===================== LOOP POR TEMPERATURAS ===========================
for expn, T_nom in zip(expns, T_nominal):
    datos = DatosProcesadosT1(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    tau, signal = datos.get_T1data(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    re = datos.espectro.real[tau.size-1]

    max_index = np.argmax(re)
    ppm_of_max = ppmAxis[max_index]
    ppmRange = [ppm_of_max - ppmIntegrationWidth / 2,
                ppm_of_max + ppmIntegrationWidth / 2]
    tau, signal = datos.get_T1data(ppmRange)

    # === PLOT ESPECTRO ===
    fig1d, ax1d = plt.subplots()
    ax1d.set_title(f"Expn {expn} — T nominal = {T_nom} °C")
    ax1d.axvspan(*ppmRange, alpha=0.2, color='b')
    ax1d.axhline(0, color='gray')
    ax1d.plot(ppmAxis, re / np.max(re), 'k', lw=2)
    ax1d.text(ppmRange[0] - 0.1, 0.8, "Regi\u00f3n de integraci\u00f3n\n(T1)", color='b')
    ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax1d.set_xlabel("Desplazamiento qu\u00edmico [ppm]")
    plt.tight_layout()
    plt.show()

    # === AJUSTE MONO ===
    tau_fit, signal_fit, residuals = datos.T1fit()
    T1_s = datos.T1params[1] / 1000
    T1mono_list.append(T1_s)

    try:
        T1_err = datos.T1stderr[1] / 1000
    except AttributeError:
        T1_err = np.std(residuals) / np.max(signal) * T1_s
    T1mono_err_list.append(T1_err)

    # === AJUSTE BI (opcional) ===
    if fit_biexponential:
        _, signal_fit_bi, _ = datos.T1fit(model='bi')
        T1short_list.append(min(datos.T1params[1], datos.T1params[3]) / 1000)
        T1long_list.append(max(datos.T1params[1], datos.T1params[3]) / 1000)
        Ashort = datos.T1params[0] if datos.T1params[1] < datos.T1params[3] else datos.T1params[2]
        Along  = datos.T1params[0] if Ashort != datos.T1params[0] else datos.T1params[2]
        P1short = Ashort / (Ashort + Along)
        P1short_list.append(P1short)

    # === SUBPLOTS ===
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(f"expn: {expn} — T nominal = {T_nom} °C", fontsize=14)

    axs[0, 0].plot(tau, signal, 'ko', label='Datos')
    axs[0, 0].plot(tau_fit, signal_fit, 'b-', label='Mono')
    if fit_biexponential:
        axs[0, 0].plot(tau_fit, signal_fit_bi, 'r--', label='Bi')
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    axs[0, 0].legend()

    axs[1, 0].plot(tau, residuals, 'ko')
    axs[1, 0].axhline(0, color='gray', linestyle='--')
    axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel='Residuos')

    axs[0, 1].plot(tau, signal, 'ko')
    axs[0, 1].plot(tau_fit, signal_fit, 'b-')
    if fit_biexponential:
        axs[0, 1].plot(tau_fit, signal_fit_bi, 'r--')
    axs[0, 1].set_xscale('log')
    axs[0, 1].set(xlabel=r'$\tau$ [ms]')
    axs[0, 1].set_yticklabels([])

    axs[1, 1].plot(tau, residuals, 'ko')
    axs[1, 1].axhline(0, color='gray', linestyle='--')
    axs[1, 1].set_xscale('log')
    axs[1, 1].set(xlabel=r'$\tau$ [ms]')
    axs[1, 1].set_yticklabels([])

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

#=======================================================================
# AJUSTE ARRHENIUS (mono)
#=======================================================================
T_real = np.array(T_real)
T1mono = np.array(T1mono_list)
T1mono_err = np.array(T1mono_err_list)

invT_K = 1 / (T_real + 273.15)
lnT1 = np.log(T1mono)
lnT1_err = T1mono_err / T1mono

def lnT1_arrhenius(invT, Ea, ln_tau0):
    R = 8.314
    return Ea / R * invT + ln_tau0

popt, pcov = curve_fit(lnT1_arrhenius, invT_K, lnT1, sigma=lnT1_err, absolute_sigma=True)
Ea_fit, ln_tau0 = popt
Ea_err = np.sqrt(np.diag(pcov))[0]
tau0 = np.exp(ln_tau0)

#%%=======================================================================
# GRAFICO ARRHENIUS FINAL
#=======================================================================
plt.figure(figsize=(7, 5))
plt.errorbar(1000*invT_K, np.exp(lnT1), yerr=lnT1_err, fmt='o', color='k', label="T1 mono")
invT_plot = np.linspace(min(invT_K)*0.95, max(invT_K)*1.05, 200)
plt.plot(1000*invT_plot, np.exp(lnT1_arrhenius(invT_plot, *popt)), 'b-', label=f'$E_a$ = {Ea_fit/1000:.2f} ± {Ea_err/1000:.2f} kJ/mol')

if fit_biexponential:
    plt.plot(1000*invT_K, T1short_list, 'rs--', label='T1 short')
    plt.plot(1000*invT_K, T1long_list, 'g^--', label='T1 long')
plt.yscale('log')
plt.xlabel(r"$1000/T$ [K$^{-1}$]")
plt.ylabel(r"$\ln(T_1)$ [s]")
plt.title("Ajuste Arrhenius: $\ln(T_1)$ vs $1000/T$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#=======================================================================
# RESUMEN
#=======================================================================
print("\nResumen:")
print("T [°C]\tT1_mono [s]\tError [s]")
for T, T1, err in zip(T_real, T1mono, T1mono_err):
    print(f"{T:.1f}\t{T1:.3f}\t{err:.3f}")

print(f"\nParámetros del ajuste Arrhenius (mono):")
print(f"Ea = {Ea_fit/1000:.2f} ± {Ea_err/1000:.2f} kJ/mol")
print(f"τ0 = {tau0:.2e} s")

if fit_biexponential:
    print("\nT1short, T1long y P1short:")
    for T, s, l, p in zip(T_real, T1short_list, T1long_list, P1short_list):
        print(f"{T:.1f} °C:\tT1short = {s:.3f} s, T1long = {l:.3f} s, P1short = {p:.2f}")

# %%
