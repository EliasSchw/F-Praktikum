import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))


e = const.e
m_eff = 1*const.m_e
epsilon_0 = const.epsilon_0
N_Dichte = 5.9*10**28
tau1 = 0.1*10**-12
tau2 = 0.01*10**-12
tau3 = 0.001*10**-12
tau4 = 0.0001*10**-12
tau5 = 3*10**-6
c = const.c
d = 500*10**-6 
k_max_semicon = 1400000
epsilon_inf = 10.3648
omega_LO = 2 * np.pi * 36.1 / const.h
omega_TO = 2 * np.pi * 33.2 / const.h
gamma = 2.5 * 100
N_DichteSemicon = 1.05 * 10**24 
m_eff_Semicon = 0.063 * const.m_e

k = np.linspace(1, k_max_semicon, 50000)
omega = k * c
epsilon_S = epsilon_inf * (1 + (omega_LO**2 - omega_TO**2) / (omega_TO**2 - omega**2 - 1j * omega * gamma))
sigma = (N_DichteSemicon * e**2 * tau5) / (m_eff_Semicon) * (1 / (1 - 1j * omega * tau5))
epsilon = epsilon_S + 1j * sigma / (omega * epsilon_0)
kappa = np.imag(np.sqrt(epsilon))
n = np.real(np.sqrt(epsilon))
beta = 4 * np.pi * kappa * omega
rh = np.abs((np.sqrt(epsilon) - 1) / (np.sqrt(epsilon) + 1))**2
R_S = rh * (1 + np.exp(-2 * beta * d) - 2 * rh * np.exp(-2 * beta * d)) / (1 - (rh**2) * np.exp(-2 * beta * d))

omega_plasma_semicon = np.sqrt((N_DichteSemicon * e**2) / (m_eff_Semicon * epsilon_0))
k_plasma_semicon = omega_plasma_semicon / c
print(k_plasma_semicon)

# Plot Reflectivity
plt.plot(k, R_S, label="Reflectivity $R_S$")
# Vertikalen Strich bei k_plasma_semicon hinzufügen
plt.axvline(x=k_plasma_semicon, color='red', linestyle='--', label="$k_{plasma}$")
plt.xlabel("k (1/cm)")
plt.ylabel("Reflectivity $R_S$")
plt.legend()
plt.title("Reflectivity vs. k")
plt.show()
