import numpy as np
import matplotlib.pyplot as plt
from Laplace import ILT
from scipy.optimize import curve_fit


def simulate_two_site_exchange(
    ppm_range=(-10, 10),
    n_points=2000,
    populations=(0.5, 0.5),       # (P1, P2)
    shifts_ppm=(-5.2, 0),         # (δ1, δ2)
    linewidths=(100, 210),        # (Δν1, Δν2) in Hz
    k_ex=0,                       # exchange rate in Hz
    nu_L=116.6e6,                  # Larmor frequency in Hz
    amp = 1.0,                # amplitude of the lineshape
    ):
    ppm = np.linspace(*ppm_range, n_points)
    freq_Hz = ppm * nu_L * 1e-6

    # Parámetros del entorno 1 (e.g., free)
    P1 = populations[0]
    delta1_Hz = shifts_ppm[0] * nu_L * 1e-6
    lw1 = linewidths[0]

    # Parámetros del entorno 2 (e.g., in-pore)
    P2 = populations[1]
    delta2_Hz = shifts_ppm[1] * nu_L * 1e-6
    lw2 = linewidths[1]

    # Funciones F1 y F2 según la ecuación S2
    F1 = P1 / (1j * (freq_Hz - delta1_Hz) - lw1 / 2 - k_ex)
    F2 = P2 / (1j * (freq_Hz - delta2_Hz) - lw2 / 2 - k_ex)

    # Línea espectral de intercambio (ecuación S1 simplificada)
    L_ex = -1j * (F1 + F2) / (1+k_ex * (F1 + F2)) / (nu_L * 1e-6)

    absorption = np.imag(L_ex) * amp
    return ppm, absorption  # Normalized intensity

# Parámetros de ejemplo (similar a los del documento, para 7Li a 7T)
ppm, spectrum = simulate_two_site_exchange(
    k_ex=100,                       # intercambio en Hz
                 # 7Li a 7T
)

plt.figure(figsize=(8, 4))
plt.plot(ppm, spectrum)
plt.gca().invert_xaxis()
plt.xlabel('Chemical shift (ppm)')
plt.ylabel('Normalized intensity')
plt.title('Two-site exchange NMR lineshape')
plt.grid(True)
plt.tight_layout()
plt.show()





ppm, spectrum_2 = simulate_two_site_exchange(
    k_ex=1,                       # intercambio en Hz
)




simulated_ppm = ppm
simulated_experiment = (spectrum + spectrum_2)

# Add Gaussian noise with SNR = 100
signal_power = np.mean(simulated_experiment ** 2)
noise_power = signal_power / 100**2
noise = np.random.normal(0, np.sqrt(noise_power), simulated_experiment.shape)
simulated_experiment += noise


plt.figure(figsize=(8, 4))
plt.plot(simulated_ppm, simulated_experiment, 'o')
plt.plot(ppm, spectrum_2, linestyle='--', color='orange', label='k_ex=1 Hz')
plt.plot(ppm, spectrum, linestyle='--', color='red', label='k_ex=100 Hz')
plt.gca().invert_xaxis()
plt.xlabel('Chemical shift (ppm)')
plt.ylabel('Normalized intensity')
plt.title('Two-site exchange NMR lineshape')
plt.grid(True)
plt.tight_layout()
plt.show()


#%% ILT==========================================================
Nilt = 50
alpha = 1e-0

# funcion de kernel
def KERNEL(ppm, k_ex):
    """
    Kernel function for the ILT.
    """
    npoints = ppm.size
    _, spectrum = simulate_two_site_exchange(
    n_points=npoints,
    k_ex=k_ex,                       # exchange rate in Hz
    )
    return spectrum


# Inicializo la clase ILT
ilt = ILT(rango=(1e-3, 1e2), kernel=KERNEL, Nilt=Nilt,
            figure=2, savepath=None)
# calculo la ILT para el conjunto de datos 1
ilt.DoTheStuff(simulated_experiment, simulated_ppm, muestra=fr"$\alpha$={ilt.alpha:.0e}")
# # calculo la ILT para el conjunto de datos 2
# ilt.DoTheStuff(ppm, spectrum, muestra="data_2")
ilt.legend()


#%%
for k in np.logspace(-1, 3, Nilt):
    plt.plot(simulated_ppm, KERNEL(simulated_ppm, k), label=f'k_ex={k:.0e} Hz', alpha=0.5)


# %%

# Define a wrapper for curve_fit that only fits k_ex
def fit_kernel(ppm, k_ex):
    return KERNEL(ppm, k_ex)

# Initial guess for k_ex
k_ex_guess = 80

# Fit the simulated data
popt, pcov = curve_fit(
    lambda ppm, k_ex: KERNEL(ppm, k_ex),
    simulated_ppm,
    simulated_experiment,
    p0=[k_ex_guess],
    bounds=(1e-1, 1e3)
)

print(f"Fitted k_ex: {popt[0]:.2f} Hz")

# Plot the fitted kernel
plt.figure(figsize=(8, 4))
plt.plot(simulated_ppm, KERNEL(simulated_ppm, popt[0]), label=f'Fitted k_ex={popt[0]:.2f} Hz', color='red')
plt.plot(simulated_ppm, simulated_experiment, 'o', label='Simulated data')
plt.gca().invert_xaxis()
plt.xlabel('Chemical shift (ppm)')
plt.ylabel('Normalized intensity')
plt.title('Fit of Two-site Exchange NMR lineshape')
plt.legend()
plt.tight_layout()
plt.show()
# %%
