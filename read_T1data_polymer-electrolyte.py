import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from scipy.optimize import curve_fit

#========================= INPUT =======================================
nexp_base = 35
T_nominal = np.arange(-20, 41, 10)
T_real = [-17.3, -9.9, -0.1, 10, 20, 30, 40]

path = r"C:\Users\Santi\OneDrive - University of Cambridge\NMRdata\500\2025-07-16_PEO-solid-electrolyte_VT/"
plotRange = [100, -100] 
# plotRange = [-25, -125]
ppmIntegrationWidth = 30

#=======================================================================
expns = [nexp_base + 10 * i for i in range(len(T_nominal))]
T1_list = []
T1_err_list = []

for expn, T_nom in zip(expns, T_nominal):
    datos = DatosProcesadosT1(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    # for === T1 ANALYSIS ===
    tau, signal = datos.get_T1data(plotRange)


    ppmAxis = datos.espectro.ppmAxis
    re = datos.espectro.real[tau.size-1]
    im = datos.espectro.imag[tau.size-1]

    max_index = np.argmax(re)
    ppm_of_max = ppmAxis[max_index]
    ppmRange = [ppm_of_max - ppmIntegrationWidth / 2,
                ppm_of_max + ppmIntegrationWidth / 2]

    # === PLOT DEL ESPECTRO ===
    r1, r2 = [np.min(ppmRange), np.max(ppmRange)]
    fig1d, ax1d = plt.subplots()
    ax1d.set_title(f"Expn {expn} — T nominal = {T_nom} °C")
    ax1d.axvspan(r1, r2, alpha=0.2, color='b')
    ax1d.axhline(0, color='gray')
    ax1d.plot(ppmAxis, re / np.max(re), 'k', lw=2)
    ax1d.text(r1 - np.abs(0.1 * r1), 0.8, "Región de integración\n(T1)", color='b')
    ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax1d.set_xlabel("Desplazamiento químico [ppm]")
    plt.tight_layout()
    plt.show()

    # === T1 ANALYSIS ===
    tau_fit, signal_fit, residuals = datos.T1fit()
    T1_s = datos.T1params[1] / 1000
    T1_list.append(T1_s)

    try:
        T1_err = datos.T1stderr[1] / 1000
    except AttributeError:
        T1_err = np.std(residuals) / np.max(signal) * T1_s
    T1_err_list.append(T1_err)

    # === Subplots ajuste T1 ===
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(f"expn: {expn} — T nominal = {T_nom} °C", fontsize=14)

    axs[0, 0].plot(tau, signal, 'ko', label='Datos')
    axs[0, 0].plot(tau_fit, signal_fit, 'b-', label='Ajuste')
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    axs[0, 0].legend()
    text = f"$T_1 =$ {datos.T1params[1]:.0f} ms\nA = {datos.T1params[0]:.2f}\n$y_0 =$ {datos.T1params[2]:.2f}"
    axs[0, 0].text(tau[-1]*0.5, (signal[-1]-signal[0])*0.15+signal[0], text)

    axs[1, 0].plot(tau, residuals, 'ko')
    axs[1, 0].axhline(0, color='gray', linestyle='--')
    axs[1, 0].set(xlabel=r'$\tau$ [ms]', ylabel='Residuos')

    axs[0, 1].plot(tau, signal, 'ko')
    axs[0, 1].plot(tau_fit, signal_fit, 'b-')
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
# AJUSTE ARRHENIUS: ln(T1) vs 1/T
#=======================================================================
T_real = np.array(T_real)
T1_list = np.array(T1_list)
T1_err_list = np.array(T1_err_list)

invT_K = 1 / (T_real + 273.15)
lnT1 = np.log(T1_list)
lnT1_err = T1_err_list / T1_list  # propagación: d(ln y) = dy / y

def lnT1_arrhenius(invT, Ea, ln_tau0):
    R = 8.314
    return Ea / R * invT + ln_tau0

popt, pcov = curve_fit(lnT1_arrhenius, invT_K, lnT1, sigma=lnT1_err, absolute_sigma=True)
Ea_fit, ln_tau0 = popt
Ea_err = np.sqrt(np.diag(pcov))[0]
tau0 = np.exp(ln_tau0)

#=======================================================================
# GRAFICO ARRHENIUS
#=======================================================================
plt.figure(figsize=(7, 5))
plt.errorbar(invT_K, lnT1, yerr=lnT1_err, fmt='o', color='k', label="Datos")
invT_plot = np.linspace(min(invT_K)*0.95, max(invT_K)*1.05, 200)
lnT1_fit = lnT1_arrhenius(invT_plot, *popt)
plt.plot(invT_plot, lnT1_fit, 'b-', label=f'Ajuste\n$E_a$ = {Ea_fit/1000:.2f} ± {Ea_err/1000:.2f} kJ/mol')

plt.xlabel(r"$1/T$ [K$^{-1}$]")
plt.ylabel(r"$\ln(T_1)$ [s]")
plt.title("Ajuste Arrhenius: $\ln(T_1)$ vs $1/T$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#=======================================================================
# RESUMEN EN TERMINAL
#=======================================================================
print("\nResumen:")
print("T [°C]\tT1 [s]\tError [s]")
for T, T1, err in zip(T_real, T1_list, T1_err_list):
    print(f"{T:.1f}\t{T1:.3f}\t{err:.3f}")

print(f"\nParámetros del ajuste Arrhenius (ln T1 vs 1/T):")
print(f"Ea = {Ea_fit/1000:.2f} ± {Ea_err/1000:.2f} kJ/mol")
print(f"τ₀ = {tau0:.2e} s")
