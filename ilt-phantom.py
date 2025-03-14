# -*- coding: utf-8 -*-
"""
Created on Wed 12/3/2025

@author: Santi
"""

import matplotlib.pyplot as plt
import numpy as np
from Datos import *
from Laplace import ILT
from scipy.optimize import curve_fit


# Define a Gaussian function
def gaussian(x, mu, sigma, amplitude):
    return amplitude * np.exp(-((np.log(x) - np.log(mu))**2) / (2 * sigma**2))

# Set SNR (adjust as needed)
SNR = 100  # Higher values → less noise, lower values → more noise

# Generate relaxation times on a logarithmic scale
relaxation_time = np.logspace(-1, 1, 100)  # From 0.1 to 10

# Use the provided parameters (mu, sigma, amplitude)
params = [
    (0.5, 0.04, 0.6),  # Gaussian centered at 0.5
    (1, 0.03, 1.0),    # Gaussian centered at 1
    (5, 0.01, 0.5)     # Gaussian centered at 5
]

# Compute the clean distribution function (sum of Gaussians)
distribution = sum(gaussian(relaxation_time, mu, sigma, amplitude) for mu, sigma, amplitude in params)

# Define a time array (linear scale)
time = np.linspace(0, 10, 300)  # From 0 to 5 with 300 points

# Compute the clean signal as a sum of exponential decays
signal = np.sum(distribution[:, np.newaxis] * np.exp(-time / relaxation_time[:, np.newaxis]), axis=0)

# Add Gaussian noise to the signal
noise = (np.max(signal) / SNR) * np.random.normal(0, 1, size=signal.shape)
signal_noisy = signal + noise  # Noisy signal

# Create the figure with two subplots
fig, ax = plt.subplots(2, 1, figsize=(8, 10))

# Plot the clean Distribution vs. Relaxation Times (scatter plot)
ax[0].scatter(relaxation_time, distribution, color="blue", s=10, label="Clean Distribution")  # Scatter plot
ax[0].set_xscale("log")  # Logarithmic x-axis
ax[0].set_xlabel("Relaxation Times", fontsize=14)
ax[0].set_ylabel("Distribution", fontsize=14)
ax[0].set_title("Distribution vs. Relaxation Times (Clean)", fontsize=16)
ax[0].legend(fontsize=12)
ax[0].grid(True, which="both", linestyle="--", linewidth=0.5)

# Plot the Noisy Signal vs. Time
ax[1].plot(time, signal_noisy, linewidth=2, color="blue", label="Noisy Signal")
ax[1].plot(time, signal, linewidth=1.5, linestyle="dashed", color="red", label="Clean Signal")  # Reference
ax[1].set_xlabel("Time", fontsize=14)
ax[1].set_ylabel("Signal", fontsize=14)
ax[1].set_title(f"Signal (Sum of Exponentials with Noise, SNR = {SNR})", fontsize=16)
ax[1].legend(fontsize=12)
ax[1].grid(True, linestyle="--", linewidth=0.5)

# Increase tick label font sizes
for axis in ax:
    axis.tick_params(axis="both", which="major", labelsize=12)

# Show the figure
plt.tight_layout()
plt.show()


#%% ILT==========================================================
Nilt = 256
alpha = 1e-2

# funcion de kernel
def T2(x, tt):
    return np.exp(-x/tt)

# Inicializo la clase ILT
ilt = ILT(rango=(1e-1, 1e1), kernel=T2, Nilt=Nilt,
            figure=2, savepath=None)
# calculo la ILT para el conjunto de datos 1
for ii in [0, -1, -2, -3,-6]:
    ilt.alpha = 10**ii
    ilt.DoTheStuff(signal_noisy, time, muestra=fr"$\alpha$={ilt.alpha:.0e}")
# # calculo la ILT para el conjunto de datos 2
# ilt.DoTheStuff(ydata2, xdata, muestra="data_2")
ilt.legend()

# %%=====================================================================
#  Fit the Noisy Signal with Exponential Models
# Exponential decay models
def exp1(t, A1, tau1):
    return A1 * np.exp(-t / tau1)

def exp2(t, A1, tau1, A2, tau2):
    return A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

def exp3(t, A1, tau1, A2, tau2, A3, tau3):
    return A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2) + A3 * np.exp(-t / tau3)
# Initial guesses
tau_guess = [0.5, 2, 5]  # Initial guesses for time constants
A_guess = [np.max(signal_noisy) / 3] * 3  # Approximate amplitudes

# Fit with one exponential
popt1, _ = curve_fit(exp1, time, signal_noisy, p0=[A_guess[0], tau_guess[0]])

# Fit with two exponentials
popt2, _ = curve_fit(exp2, time, signal_noisy, p0=[A_guess[0], tau_guess[0], A_guess[1], tau_guess[1]])

# Fit with three exponentials
popt3, _ = curve_fit(exp3, time, signal_noisy, p0=[A_guess[0], tau_guess[0], A_guess[1], tau_guess[1], A_guess[2], tau_guess[2]])

# Compute fitted signals
signal_fit1 = exp1(time, *popt1)
signal_fit2 = exp2(time, *popt2)
signal_fit3 = exp3(time, *popt3)

# Compute residuals
residual1 = signal_noisy - signal_fit1
residual2 = signal_noisy - signal_fit2
residual3 = signal_noisy - signal_fit3


fig, ax = plt.subplots(3, 1, figsize=(8, 12))

# Plot the clean Distribution vs. Relaxation Times
ax[0].scatter(relaxation_time, distribution, color="blue", s=10, label="Clean Distribution")
ax[0].set_xscale("log")  # Logarithmic x-axis
ax[0].set_xlabel("Relaxation Times", fontsize=14)
ax[0].set_ylabel("Distribution", fontsize=14)
ax[0].set_title("Distribution vs. Relaxation Times", fontsize=16)
ax[0].legend(fontsize=12)
ax[0].grid(True, which="both", linestyle="--", linewidth=0.5)

# Plot the Noisy Signal and Fits
ax[1].plot(time, signal_noisy, "o", markersize=3, color="blue", alpha=0.6, label="Noisy Signal")
ax[1].plot(time, signal_fit1, "--", color="green", label="Fit: 1 Exp")
ax[1].plot(time, signal_fit2, "--", color="orange", label="Fit: 2 Exp")
ax[1].plot(time, signal_fit3, "--", color="red", label="Fit: 3 Exp")
ax[1].set_xlabel("Time", fontsize=14)
ax[1].set_ylabel("Signal", fontsize=14)
ax[1].set_title("Fitting of Noisy Signal", fontsize=16)
ax[1].legend(fontsize=12)
ax[1].grid(True, linestyle="--", linewidth=0.5)

# Plot the Residuals
ax[2].plot(time, residual1, "-", color="green", label="Residual: 1 Exp")
ax[2].plot(time, residual2, "-", color="orange", label="Residual: 2 Exp")
ax[2].plot(time, residual3, "-", color="red", label="Residual: 3 Exp")
ax[2].axhline(0, color="black", linestyle="--", linewidth=1)
ax[2].set_xlabel("Time", fontsize=14)
ax[2].set_ylabel("Residuals", fontsize=14)
ax[2].set_title("Residuals of Fits", fontsize=16)
ax[2].legend(fontsize=12)
ax[2].grid(True, linestyle="--", linewidth=0.5)

# Increase tick label font sizes
for axis in ax:
    axis.tick_params(axis="both", which="major", labelsize=12)

# Show the figure
plt.tight_layout()
plt.show()

# %%
