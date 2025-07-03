import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv  # Bessel function of the first kind



def redor_curve(r, spin_speed, N):
    """
    Calculate the recoupling curve for a given internuclear distance and spinning speed.
    
    Parameters:
    r (float): Internuclear distance in nm.
    spin_speed (float): Spinning speed in Hz.
    N (array-like): Array of rotor cycles.
    
    Returns:
    tuple: Recoupling times and normalized signal intensities.
    """
    # Convert internuclear distance to meters
    R = r * 1e-9  # nm to m
    mu = 4 * np.pi * 1e-7  # magnetic constant (H/m)
    gamma_7Li = 103.962e6  # rad/(s*T)
    gamma_F19 = 251.662e6  # rad/(s*T)
    h_bar = 1.05457e-34  # J*s

    # Calculate dipolar coupling constant
    D = mu * gamma_7Li * gamma_F19 * h_bar / (8 * np.pi**2 * R**3)

    # Calculate recoupling time
    Tr = 1 / spin_speed  # rotor period in seconds
    NTr = Tr * N  # recoupling time
    _lambda = N * Tr * D

    # Calculate Bessel function
    S_S0 = np.sqrt(1/8) * np.pi * jv(0.25, np.sqrt(2) * _lambda) * jv(-0.25, np.sqrt(2) * _lambda)

    return NTr, S_S0


if __name__ == "__main__":
    plt.figure(1)
    for r in np.arange(0.25, 1.51, 0.25):  # Example range of internuclear distances
        
        # r:  internuclear distance in nm
        spin_speed = 14000  # spinning speed in Hz
        N = np.arange(1, 65)  # rotor cycles from 1 to 100

        NTr, S_S0 = redor_curve(r, spin_speed, N)

        # ====== Plot ===========
        plt.plot(NTr*1000, S_S0, label=f"r = {r} nm")
    plt.xlabel('Recoupling Time (ms)')
    plt.ylabel('S/S₀')
    plt.title('Recoupling Curve')
    plt.grid(True)
    plt.legend()
    plt.show()
