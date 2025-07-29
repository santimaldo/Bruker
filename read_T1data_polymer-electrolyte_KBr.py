import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from scipy.interpolate import interp1d

#========================= INPUT =======================================
nexp_base = 30
plotRange = [100, -100]

T_nominal = np.arange(-20, 41, 10)
T_real = [-17.8, -9.9, 0, 10, 20, 30, 40]

path = r"C:\\Users\\Santi\\OneDrive - University of Cambridge\\NMRdata\\500\\2025-07-16_PEO-solid-electrolyte_VT/"
ppmIntegrationWidth = 30

#========================= INIT ========================================
expns_base = [nexp_base + 10 * i for i in range(len(T_nominal))]
expns_final = [e + 7 for e in expns_base]

T1mono_base_list, T1mono_base_err_list = [], []
T1mono_final_list, T1mono_final_err_list = [], []

ppm_max_base_list = []
ppm_max_final_list = []

#===================== FUNCIONES AUXILIARES ============================
def procesar_exp(expn, T_nom):
    datos = DatosProcesadosT1(f'{path}/{expn}/')
    datos.espectro.ppmSelect(plotRange)

    tau, signal = datos.get_T1data(plotRange)

    ppmAxis = datos.espectro.ppmAxis
    re = datos.espectro.real[-1]

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
    ax1d.text(ppmRange[0] - 0.1, 0.8, "Región de integración\n(T1)", color='b')
    ax1d.set_xlim(np.max(ppmAxis), np.min(ppmAxis))
    ax1d.set_xlabel("Desplazamiento químico [ppm]")
    plt.tight_layout()
    plt.show()

    # === AJUSTE MONO ===
    tau_fit, signal_fit, residuals = datos.T1fit()
    T1_s = datos.T1params[1] / 1000

    try:
        T1_err = datos.T1stderr[1] / 1000
    except AttributeError:
        T1_err = np.std(residuals) / np.max(signal) * T1_s

    # === PLOT AJUSTE ===
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(f"Expn: {expn} — T nominal = {T_nom} °C", fontsize=14)

    axs[0, 0].plot(tau, signal, 'ko', label='Datos')
    axs[0, 0].plot(tau_fit, signal_fit, 'b-', label='Mono')
    axs[0, 0].set(xlabel=r'$\tau$ [ms]', ylabel=r'$S_{norm}$')
    axs[0, 0].legend()

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

    return T1_s, T1_err, ppm_of_max

#===================== LOOP POR PARES DE EXPERIMENTOS =================
for T_nom, exp_base, exp_final in zip(T_nominal, expns_base, expns_final):
    T1_b, err_b, ppm_b = procesar_exp(exp_base, T_nom)
    T1_f, err_f, ppm_f = procesar_exp(exp_final, T_nom)

    T1mono_base_list.append(T1_b)
    T1mono_base_err_list.append(err_b)
    ppm_max_base_list.append(ppm_b)

    T1mono_final_list.append(T1_f)
    T1mono_final_err_list.append(err_f)
    ppm_max_final_list.append(ppm_f)

#==================== GRAFICO T1 VS TEMPERATURA =========================
plt.figure(figsize=(7, 5))
plt.errorbar(T_real, T1mono_base_list, yerr=T1mono_base_err_list, fmt='o-', color='k', label='Inicial (expn base)')
plt.errorbar(T_real, T1mono_final_list, yerr=T1mono_final_err_list, fmt='s--', color='r', label='Final (expn base + 7)')
plt.xlabel("Temperatura [°C]")
plt.ylabel(r"$T_1$ [s]")
plt.title("Relajación T1 vs Temperatura (ajuste monoexponencial)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#==================== GRAFICO DESPLAZAMIENTO QUIMICO ====================
plt.figure(figsize=(7, 5))
plt.plot(T_real, ppm_max_base_list, 'o-', color='k', label='Inicial (expn base)')
plt.plot(T_real, ppm_max_final_list, 's--', color='r', label='Final (expn base + 7)')
plt.xlabel("Temperatura [°C]")
plt.ylabel("Desplazamiento químico [ppm]")
plt.title("Desplazamiento químico del máximo del espectro")
plt.gca().invert_xaxis()
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#==================== RESUMEN ============================================
print("\nResumen T1 y ppm máximo:")
print("T [°C]\tT1_ini [s]\tErr_ini [s]\tppm_ini\tT1_fin [s]\tErr_fin [s]\tppm_fin")
for T, Tb, eb, pb, Tf, ef, pf in zip(T_real,
                                      T1mono_base_list, T1mono_base_err_list, ppm_max_base_list,
                                      T1mono_final_list, T1mono_final_err_list, ppm_max_final_list):
    print(f"{T:.1f}\t{Tb:.3f}\t{eb:.3f}\t{pb:.2f}\t{Tf:.3f}\t{ef:.3f}\t{pf:.2f}")

#%%=====================================
#======================================= Temperature calibration


############# IMPLEMENTAR

# # --------------------------------------
# # Thurber empirical model for T1 vs T
# # --------------------------------------
# def thurber_model(T):
#     """
#     Empirical model by Thurber for T1 relaxation time as a function of temperature (K).
#     """
#     return 0.0145 + 5330 * T**-2 + 1.42e7 * T**-4 + 2.48e9 * T**-6


# T1_list = T1mono_base_list

# # Inverse interpolation: T1 → T
# T = np.linspace(290, 350, 1000)
# T1_model = thurber_model(T)
# T1_to_T = interp1d(T1_model, T, kind='linear', fill_value='extrapolate')

# # Temperature estimated from T1
# T_from_T1_C = T1_to_T(T1_list / 1000) - 273.15  # in °C

# plt.figure(figsize=(6, 4))
# plt.plot(time_minutes, T_from_T1_C, 'o-', label='Estimated T from $T_1$')
# plt.xlabel("Time [min]")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Temperature estimated from $T_1$")
# plt.grid(True)
# plt.tight_layout()

# # Plot: T1 vs Temperature (model vs experimental)
# fig, ax = plt.subplots(figsize=(8, 4))
# T_plot = np.linspace(20, 296, 1000)  # Temperature range fixed as requested
# T1_plot = thurber_model(T_plot)

# ax.plot(T_plot, T1_plot, label='Thurber model')
# ax.plot(T1_to_T(T1_list / 1000), T1_list / 1000, 'o-', label='Experimental data')
# ax.set_xlabel("Temperature [K]")
# ax.set_ylabel(r"$T_1$ [s]")
# ax.set_yscale('log')
# ax.set_title(r"$T_1$ vs Temperature")
# ax.legend()
# ax.grid(True, which='both', linestyle='--', linewidth=0.5)
# fig.tight_layout()


# # --- Chemical shift to temperature calibration ---
# ppm_valid = ppm_max_base_list

# T_chemshift_slope = -0.025  # ppm/K
# T0 = T1_to_T(T1_list[0] / 1000)
# ppm0 = ppm_valid[0]
# offset_ppm = ppm0 - T_chemshift_slope * T0

# def chemshift_to_T(ppm):
#     """
#     Convert chemical shift (ppm) to temperature (K).
#     """
#     return (ppm - offset_ppm) / T_chemshift_slope

# # Estimated temperature from chemical shift
# T_from_shift_C = chemshift_to_T(ppm_valid) - 273.15

# plt.figure(figsize=(6, 4))
# plt.plot(T_real, T_from_shift_C, 'o-', label='Estimated T from chemical shift')
# plt.xlabel("T")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Temperature estimated from chemical shift")
# plt.grid(True)
# plt.tight_layout()

# # --- Combined plot: temperature estimates from both methods ---
# plt.figure(figsize=(7, 4))
# plt.plot(T_real, T_from_T1_C, 'o-', label='From $T_1$')
# plt.plot(T_real, T_from_shift_C, 'o-', label='From chemical shift')
# plt.xlabel("T")
# plt.ylabel("Estimated temperature [°C]")
# plt.title("Comparison of temperature estimation methods")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# # %%
