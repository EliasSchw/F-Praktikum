import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))


e = const.e
m_eff = 1*const.m_e
epsilon_0 = 1.66818
N_Dichte = 5.9*10**28
tau1 = 0.1
tau2 = 0.01
tau3 = 0.001
tau4 = 0.0001
c=const.c
d = 500*10**-5 #dicke si undoped
k_max = 1*10**8
k_max_semicon = 1500
epsilon_r = 10.3648
omega_LO = 292
omega_TO = 268
gamma = 2.5
N_DichteSemicon = 1.05
m_eff_Semicon = 0.067
esquared = 28202.2


def calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau):
    omega = k  
    epsilon_S = epsilon_inf * (1 + (omega_LO**2 - omega_TO**2) / (omega_TO**2 - omega**2 - 1j * omega * gamma))
    sigma = (N_DichteSemicon * esquared * tau) / (m_eff_Semicon) * (1 / (1 - 1j * omega * tau))
    epsilon = epsilon_S + 1j * sigma / (omega * epsilon_0)
    kappa = np.imag(np.sqrt(epsilon))
    n = np.real(np.sqrt(epsilon))
    beta = 4 * np.pi * kappa * omega
    rh = np.abs((np.sqrt(epsilon) - 1) / (np.sqrt(epsilon) + 1))**2
    R_S = rh * (1 + np.exp(-2 * beta * d) - 2 * rh * np.exp(-2 * beta * d)) / (1 - (rh**2) * np.exp(-2 * beta * d))
    return R_S

def calculateKPlasmaSemiconductor(N_DichteSemicon, esquared, m_eff_Semicon, epsilon_0, epsilon_r):
    omega_plasma_semicon = np.sqrt(N_DichteSemicon * esquared / (m_eff_Semicon * epsilon_0 * epsilon_r))
    k_plasma_semicon = omega_plasma_semicon  # Keine zusätzliche Umrechnung
    return k_plasma_semicon

def plotReflectivitySemiconductor(d, tau1, tau2, tau3, tau4, k_max_semicon):
    k = np.linspace(1, k_max_semicon, 50000)  # k wird hier definiert
    k_plasma_semicon = calculateKPlasmaSemiconductor(N_DichteSemicon, esquared, m_eff_Semicon, epsilon_0, epsilon_r)

    taus = [tau1, tau2, tau3, tau4]
    fig, ax = plt.subplots()
    tausGray = np.logspace(np.log10(tau1), np.log10(tau4), 25)
    for tauGray in tausGray:
        R_S_gray = calculateReflectivitySemiconductor(epsilon_r, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tauGray)
        ax.plot(k / 1000, R_S_gray, color='gray', linewidth=0.5, alpha=0.7)  # Graue Linien

    for tau in taus:
        R_S = calculateReflectivitySemiconductor(epsilon_r, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau)
        ax.plot(k / 1000, R_S, label=f'τ={tau} ps')  # Hauptkurven

    # Vertikale gestrichelte Linie bei k_plasma_semicon
    plt.axvline(x=k_plasma_semicon / 1000, color='black', linestyle='--', label=r'$k_{plasma}$')

    # Achsentitel und Bereich anpassen
    plt.xlabel(r'wave number k [$10^3 \, \mathrm{m}^{-1}$]')
    plt.ylabel('reflectivity R')
    plt.title('Theoretical reflectivity of a semiconductor')
    plt.legend()
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.xlim(left=0, right=1.5)
    plt.ylim(bottom=0, top=1)

    plt.savefig('Paper/Images/semi.png', dpi=400)
    from PIL import Image
    Image.open("Paper/Images/semi.png").show()
    plt.clf()


plotReflectivitySemiconductor(d, tau1, tau2, tau3, tau4, k_max_semicon)



