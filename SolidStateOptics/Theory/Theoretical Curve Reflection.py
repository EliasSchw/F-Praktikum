import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

# Konstanten
e = const.e
m_eff = 1*const.m_e
epsilon_0 = const.epsilon_0
N_Dichte = 5.9*10**28
tau = 0.1*10**-12
c=const.c
d = 100 *10**-5 #dicke si undoped
k_max = 1*10**8



k = np.linspace(0.001,k_max,1000)
omega = c*k
omega_p = np.sqrt(N_Dichte*e**2/(m_eff*epsilon_0))
epsilon = 1-omega_p**2*tau**2/(1+omega**2*tau**2)+1j*omega_p**2*tau/omega/(1+omega**2*tau**2)
N = np.sqrt(epsilon)
#R_Halbraum = np.abs((N - 1)/(N + 1))**2


#Gleichung 1.60 skript
expo = 1j*2*omega*N*d/c
r=(np.exp(expo)-1)*(1-N)/(np.exp(expo)*(1-N)-(1+N))
R_FabryPerot = np.abs(r)**2


fig, ax = plt.subplots()
        
ax.plot(k,R_FabryPerot)

#ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r'$\Gamma$', 'X', 'M', 'Y', r'$\Gamma$'], size=14)

ax.set_ylabel('Reflectivity R')
ax.set_xlabel('wavenumber k /' + r'$m^{-1}$')
ax.set_title('Reflectivity of a Fabry-Perot-Interferometer')
#ax.set_xlim(left=0, right=2)
#ax.set_ylim(bottom=-3, top=3)
ax.grid()

plt.savefig('Paper/Images/foo.png', dpi=400)
from PIL import Image
Image.open("Paper/Images/foo.png").show()
plt.clf()



# plt.plot(k,R_FabryPerot)
# plt.show()
