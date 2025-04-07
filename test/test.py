import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from macroswriter import writeLatexMacro

e = const.e
m_eff = 1 * const.m_e
epsilon_0 = const.epsilon_0
tau2 = 0.1 * 10**-12
tau3 = 0.01 * 10**-12
tau4 = 0.001 * 10**-12
c = const.c
d = 500 * 10**-6
k_max_semicon = 1500000
epsilon_inf = 10.89
omega_LO = 2 * np.pi * 36.1 * 10**-3 / const.h
omega_TO = 2 * np.pi * 33.2 * 10**-3 / const.h
gamma = 2.5 * 100
N_DichteSemicon = 1.05 * 10**24
m_eff_Semicon = 0.063 * const.m_e

k = np.linspace(1, k_max_semicon, 50000)

def calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau):
    nu = 2*np.pi*k
    omega = 2*np.pi * c
    epsilon_S = epsilon_inf * (1 + (omega_LO**2 - omega_TO**2) / (omega_TO**2 - omega**2 - 1j * omega * gamma))
    sigma = (N_DichteSemicon * e**2 * tau) / (m_eff_Semicon) * (1 / (1 - 1j * omega * tau))
    epsilon = epsilon_S + 1j * sigma / (omega * epsilon_0)
    kappa = np.imag(np.sqrt(epsilon))
    n = np.real(np.sqrt(epsilon))
    beta = 4 * np.pi * kappa * omega
    rh = np.abs((np.sqrt(epsilon) - 1) / (np.sqrt(epsilon) + 1))**2
    R_S = rh + ((1-rh)**2*rh*np.exp(-2*beta*d))/(1-rh**1*np.exp(-2*beta*d))
    return R_S

def calculateKPlasmaSemiconductor(N_DichteSemicon, m_eff_Semicon, epsilon_0, epsilon_inf):
    omega_plasma_semicon = np.sqrt((N_DichteSemicon * e**2) / (m_eff_Semicon * epsilon_0 * epsilon_inf))
    k_plasma_semicon = omega_plasma_semicon / c
    return k_plasma_semicon

def plotReflectivitySemiconductor(d, tau2, tau3, tau4, k_max_semicon):
    k = np.linspace(0.1, k_max_semicon, 500000)  # k wird hier definiert
    k_plasma_semicon = calculateKPlasmaSemiconductor(N_DichteSemicon, m_eff_Semicon, epsilon_0, epsilon_inf)

    taus = [tau2, tau3, tau4]
    fig, ax = plt.subplots()
    tausGray = np.logspace(np.log10(tau2), np.log10(tau4), 25)
    for tauGray in tausGray:
        R_S_gray = calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tauGray)
        ax.plot(k, R_S_gray, color='gray', linewidth=0.5, alpha=0.7)  # Graue Linien

    for tau in taus:
        R_S = calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau)
        ax.plot(k, R_S, label=f'τ={tau*10**12} ps')  # Hauptkurven

    # Vertikale gestrichelte Linie bei k_plasma_semicon
    plt.axvline(x=k_plasma_semicon, color='black', linestyle='--', label=r'$k_{plasma}$')

    # Achsentitel und Bereich anpassen
    plt.xlabel(r'wave number k [$10^3 \, \mathrm{m}^{-1}$]')
    plt.ylabel('reflectivity R')
    plt.title('Theoretical reflectivity of a semiconductor')
    plt.legend()
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.xlim(left=0, right=1.5 * 10**6)
    plt.ylim(bottom=0, top=1)

    plt.savefig('Paper/Images/semi.png', dpi=400)
    from PIL import Image
    Image.open("Paper/Images/semi.png").show()
    plt.clf()

plotReflectivitySemiconductor(d, tau2, tau3, tau4, k_max_semicon)


