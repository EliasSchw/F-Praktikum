import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

# Konstanten
e = const.e
m_eff = 1*const.m_e
epsilon_0 = const.epsilon_0
N_Dichte = 5.9*10**28
tau = 0.1*10**-12
k = np.linspace(0.1,140*1000*100*10,1000)
c=const.c
d = 530 *10**-6 #dicke si undoped


#Elias
omega = c*k
omega_p = np.sqrt(N_Dichte*e**2/(m_eff*epsilon_0))
epsilon = 1-omega_p**2*tau**2/(1+omega**2*tau**2)+1j*omega_p**2*tau/omega/(1+omega**2*tau**2)
N = np.sqrt(epsilon)
#R_Halbraum = np.abs((N - 1)/(N + 1))**2


#Gleichung 1.60 skript
expo = 1j*2*omega*N*d/c
r=(np.exp(expo)-1)*(1-N)/(np.exp(expo)*(1-N)-(1+N))
R_FabryPerot = np.abs(r)**2


plt.plot(k,R_FabryPerot)
plt.show()
