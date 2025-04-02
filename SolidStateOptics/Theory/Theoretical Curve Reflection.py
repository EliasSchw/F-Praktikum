import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt


import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from macroswriter import writeLatexMacro  # Import macroswriter module

# Konstanten
e = const.e
m_eff = 1*const.m_e
epsilon_0 = const.epsilon_0
N_Dichte = 5.9*10**28
tau1 = 0.1*10**-12
tau2 = 0.01*10**-12
tau3 = 0.001*10**-12
tau4 = 0.0001*10**-12
tau5 = 0.00001*10**-12
c=const.c
d = 500*10**-5 #dicke si undoped
k_max = 1*10**8
k_max_semicon = 1400000
epsilon_inf = 10.3648
omega_LO = 292*100 
omega_TO = 268*100 
gamma = 2.5*100
N_DichteSemicon = 1.05*10**24 
m_eff_Semicon = 0.067*const.m_e



def calculateReflectivity(d, tau, k):
    omega = c*k
    omega_p = np.sqrt(N_Dichte*e**2/(m_eff*epsilon_0))
    epsilon = 1-omega_p**2*tau**2/(1+omega**2*tau**2)+1j*omega_p**2*tau/omega/(1+omega**2*tau**2)
    N = np.sqrt(epsilon)

    #Gleichung 1.60 skript
    expo = 1j*2*omega*N*d/c
    r=(np.exp(expo)-1)*(1-N)/(np.exp(expo)*(1-N)-(1+N))
    R_FabryPerot = np.abs(r)**2
    return R_FabryPerot

# Berechnung der Plasmafrequenz und des zugehörigen Wellenvektors für den ersten Plot
omega_plasma = np.sqrt(N_Dichte * e**2 / (m_eff * epsilon_0))
k_plasma = omega_plasma / c
def plotReflectivity(d, tau1, tau2, tau3, tau4, tau5, k_max):
    fig, ax = plt.subplots()
    
    k = np.linspace(0, k_max, 10000)

    taus = np.logspace(np.log10(tau1), np.log10(tau5), 50)  # Generate 50 tau values logarithmically spaced
    for tau in taus:
        ax.plot(k, calculateReflectivity(d, tau, k), color='gray', linewidth=0.5, alpha=0.7)
        
    ax.plot(k, calculateReflectivity(d, tau5, k), linewidth=3)
    ax.plot(k, calculateReflectivity(d, tau4, k), linewidth=3)
    ax.plot(k, calculateReflectivity(d, tau3, k), linewidth=3)
    ax.plot(k, calculateReflectivity(d, tau2, k), linewidth=3)
    ax.plot(k, calculateReflectivity(d, tau1, k), linewidth=3)
    #ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r'$\Gamma$', 'X', 'M', 'Y', r'$\Gamma$'], size=14)

    ax.set_ylabel('Reflectivity R')
    ax.set_xlabel('wavenumber k / ' + r'$m^{-1}$')
    ax.set_title('Reflectivity of a Fabry-Perot-Interferometer')
    ax.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    ax.set_xlim(left=-0.01*k_max, right=k_max)
    ax.set_ylim(bottom=-0.05, top=1.05)
    #ax.grid()

    plt.savefig('Paper/Images/foo.png', dpi=400)
    from PIL import Image
    #Image.open("Paper/Images/foo.png").show()
    plt.clf()
    plt.close()

plotReflectivity(d, tau1, tau2, tau3, tau4, tau5, k_max)


writeLatexMacro('N_Dichte',np.round(N_Dichte,4))

#Lukas

def calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau):
    omega = c * k
    epsilon_S = epsilon_inf * (1 + (omega_LO**2 - omega_TO**2) / (omega_TO**2 - omega**2 - 1j * omega * gamma))
    sigma = (N_DichteSemicon * e**2 * tau) / (m_eff_Semicon) * (1 / (1 - 1j * omega * tau))
    epsilon = epsilon_S + 1j * sigma / (omega * epsilon_0)
    N_S = np.sqrt(epsilon)
    expo = 1j * 2 * omega * N_S * d / c
    r = (np.exp(expo) - 1) * (1 - N_S) / (np.exp(expo) * (1 - N_S) - (1 + N_S))
    R_S = np.abs(r)**2  # Normierung entfernt
    return R_S

# Berechnung der Plasmafrequenz und des zugehörigen Wellenvektors für den Halbleiter-Plot
omega_plasma_semicon = np.sqrt(N_DichteSemicon * e**2 / (m_eff_Semicon * epsilon_0))
k_plasma_semicon = omega_plasma_semicon / c

# Entfernen der Plasmafrequenzberechnung und der vertikalen Linie aus dem zweiten Plot
k = np.linspace(1, k_max_semicon, 50000)  # Startwert auf 0.1 gesetzt, höhere Auflösung

taus = [tau1, tau2, tau3, tau4]
for tau in taus:
    R_S = calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau)
    # Überprüfung, ob R_S gültige Werte enthält
    if np.all(np.isnan(R_S)) or np.all(R_S == 0):
        print(f"Warnung: Alle Werte von R_S für τ={tau} sind ungültig oder Null.")
    else:
        plt.plot(k, R_S, label=f'τ={tau*1e12:.4f} ps')

# Vertikale gestrichelte Linie bei k_plasma_semicon
plt.axvline(x=k_plasma_semicon, color='black', linestyle='--', label=r'$k_{plasma}$')

# Achsentitel, Legende und Speichern
plt.xlabel('wave number k / ' + r'$m^{-1}$')
plt.ylabel('reflectivity R')
plt.title('Theoretical reflectivity of a semiconductor')
plt.legend()
plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
plt.xlim(left=0, right=k_max_semicon)
plt.ylim(bottom=0, top=1)

plt.savefig('Paper/Images/semi.png', dpi=400)
from PIL import Image
Image.open("Paper/Images/semi.png").show()
plt.clf()


writeLatexMacro('N_Dichte',np.round(N_Dichte,4))

#Lukas

def calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau):
    omega = c * k
    epsilon_S = epsilon_inf * (1 + (omega_LO**2 - omega_TO**2) / (omega_TO**2 - omega**2 - 1j * omega * gamma))
    sigma = (N_DichteSemicon * e**2 * tau) / (m_eff_Semicon) * (1 / (1 - 1j * omega * tau))
    epsilon = epsilon_S + 1j * sigma / (omega * epsilon_0)
    N_S = np.sqrt(epsilon)
    expo = 1j * 2 * omega * N_S * d / c
    r = (np.exp(expo) - 1) * (1 - N_S) / (np.exp(expo) * (1 - N_S) - (1 + N_S))
    R_S = np.abs(r)**2  # Normierung entfernt
    return R_S

# Berechnung der Plasmafrequenz und des zugehörigen Wellenvektors für den Halbleiter-Plot
omega_plasma_semicon = np.sqrt(N_DichteSemicon * e**2 / (m_eff_Semicon * epsilon_0))
k_plasma_semicon = omega_plasma_semicon / c

# Entfernen der Plasmafrequenzberechnung und der vertikalen Linie aus dem zweiten Plot
k = np.linspace(1, k_max_semicon, 50000)  # Startwert auf 0.1 gesetzt, höhere Auflösung

taus = [tau1, tau2, tau3, tau4]
for tau in taus:
    R_S = calculateReflectivitySemiconductor(epsilon_inf, omega_LO, omega_TO, gamma, d, k, N_DichteSemicon, tau)
    # Überprüfung, ob R_S gültige Werte enthält
    if np.all(np.isnan(R_S)) or np.all(R_S == 0):
        print(f"Warnung: Alle Werte von R_S für τ={tau} sind ungültig oder Null.")
    else:
        plt.plot(k, R_S, label=f'τ={tau*1e12:.4f} ps')

# Vertikale gestrichelte Linie bei k_plasma_semicon
plt.axvline(x=k_plasma_semicon, color='black', linestyle='--', label=r'$k_{plasma}$')

# Achsentitel, Legende und Speichern
plt.xlabel('wave number k / ' + r'$m^{-1}$')
plt.ylabel('reflectivity R')
plt.title('Theoretical reflectivity of a semiconductor')
plt.legend()
plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
plt.xlim(left=0, right=k_max_semicon)
plt.ylim(bottom=0, top=1)

plt.savefig('Paper/Images/semi.png', dpi=400)
from PIL import Image
Image.open("Paper/Images/semi.png").show()
plt.clf()
