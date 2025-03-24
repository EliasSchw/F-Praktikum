import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

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
d = 100 *10**-5 #dicke si undoped
k_max = 1*10**8





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


fig, ax = plt.subplots()
        
k = np.linspace(0, k_max, 10000)

taus = np.logspace(np.log10(tau1), np.log10(tau5), 50)  # Generate 50 tau values logarithmically spaced
for tau in taus:
    ax.plot(k, calculateReflectivity(d, tau, k), color='gray', linewidth=0.5, alpha=0.7)
    
ax.plot(k, calculateReflectivity(d, tau1, k), linewidth=3)
ax.plot(k, calculateReflectivity(d, tau2, k), linewidth=3)
ax.plot(k, calculateReflectivity(d, tau3, k), linewidth=3)
ax.plot(k, calculateReflectivity(d, tau4, k), linewidth=3)
ax.plot(k, calculateReflectivity(d, tau5, k), linewidth=3)





#ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r'$\Gamma$', 'X', 'M', 'Y', r'$\Gamma$'], size=14)

ax.set_ylabel('Reflectivity R')
ax.set_xlabel('wavenumber k / ' + r'$m^{-1}$')
ax.set_title('Reflectivity of a Fabry-Perot-Interferometer')
ax.tick_params(axis='both', direction='in', which='both', top=True, right=True)
ax.set_xlim(left=0, right=k_max)
ax.set_ylim(bottom=0, top=1)
#ax.grid()

plt.savefig('Paper/Images/foo.png', dpi=400)
from PIL import Image
Image.open("Paper/Images/foo.png").show()
plt.clf()



# plt.plot(k,R_FabryPerot)
# plt.show()
