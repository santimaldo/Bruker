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
    #print(f"Internuclear distance: {r} nm, Dipolar coupling D: {D:.2e} Hz")
    # Calculate recoupling time
    Tr = 1 / spin_speed  # rotor period in seconds
    NTr = Tr * N  # recoupling time
    _lambda = N * Tr * D

    # Calculate Bessel function
    S_S0 = np.sqrt(1/8) * np.pi * jv(0.25, np.sqrt(2) * _lambda) * jv(-0.25, np.sqrt(2) * _lambda)

    #relative difference
    # rd = 1 - (jv(0, np.sqrt(2)*_lambda))**2
    # for k in range(1, 11):
    #     rd += 2/(16*k**2 - 1) * (jv(k, np.sqrt(2)*_lambda))**2

    # S_S0 = 1 - rd
    return NTr, S_S0


def DeltaS_quadrupolar(r, spin_speed, N, Nterms=10):

    # Convert internuclear distance to meters
    R = r * 1e-9  # nm to m
    mu = 4 * np.pi * 1e-7  # magnetic constant (H/m)
    gamma_7Li = 103.962e6  # rad/(s*T)
    gamma_F19 = 251.662e6  # rad/(s*T)
    h_bar = 1.05457e-34  # J*s

    # Calculate dipolar coupling constant
    D = mu * gamma_7Li * gamma_F19 * h_bar / (8 * np.pi**2 * R**3)
    #print(f"Internuclear distance: {r} nm, Dipolar coupling D: {D:.2e} Hz")
    # Calculate recoupling time
    Tr = 1 / spin_speed  # rotor period in seconds
    NTr = Tr * N  # recoupling time
    _lambda = N * Tr * D

    I = 1.5 # Spin quantum number for 7Li
    rd = 0
    for m in np.arange(-I, I+1):
        rd += 1 - (jv(0, 2*np.sqrt(2)*np.abs(m)*_lambda))**2
        for k in range(1, Nterms+1):
            rd += 2/(16*k**2 - 1) * (jv(k, 2*np.sqrt(2)*np.abs(m)*_lambda))**2
    rd = rd/(2*I + 1)


    S_S0 = 1 - rd
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


    plt.figure(2)
    colors = plt.cm.tab10(np.linspace(0, 1, 6))  # Create color palette for 6 unique colors
    
    for i, r in enumerate(np.arange(0.25, 1.51, 0.25)):  
        spin_speed = 14000  
        N = np.arange(1, 65)  

        NTr, S_S0 = redor_curve(r, spin_speed, N)
        plt.plot(NTr*1000, S_S0, color=colors[i], label=f"r = {r} nm (dipolar)")
        
        NTr_q, S_S0_q = DeltaS_quadrupolar(r, spin_speed, N)
        plt.plot(NTr_q*1000, S_S0_q, '--', color=colors[i], label=f"r = {r} nm (quad)")
    
    plt.xlabel('Recoupling Time (ms)')
    plt.ylabel('S/S₀')
    plt.title('Recoupling Curve')
    plt.grid(True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()
