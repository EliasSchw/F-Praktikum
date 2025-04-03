from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
import numpy as np
import scipy.constants as const
from scipy.optimize import fsolve
from Glätteisen import bügeln

e = const.e  # Elementary charge in Coulombs
hbar = const.hbar  # Reduced Planck's constant in J.s
m = const.m_e  # Electron mass in kg
epsilon_0 = const.epsilon_0  # Vacuum permittivity in F/m

#moch willkürliche Werte
mu_GaAs = 0.036872 * const.m_e  # Effective mass of electron in GaAs (müssen quelle finden!!)
m_e = const.m_e  # Electron mass in kg
factor = e*e*(2*mu_GaAs)**(3/2)*2*np.pi/(epsilon_0*m_e**2*hbar**3)    
d = 0.1e-6  # Thickness of the sample in m



#calculate n

initial_guess = [1,1] # Initial guess for beta and R

def T_fabry_perot(vars, T_fabry_perot_value):
    beta, R = vars
    return T_fabry_perot_value - (1-R)**2*np.exp(-beta*d)/(1-R**2*np.exp(-2*beta*d))

def R_fabry_perot(vars, R_fabry_perot_value):
    beta, R = vars
    return R_fabry_perot_value -(R + (1-R)**2*R*np.exp(-2*beta*d)/(1-R**2*np.exp(-2*beta*d)))

def equations(vars, T_fabry_perot_value, R_fabry_perot_value):
    return [T_fabry_perot(vars, T_fabry_perot_value), R_fabry_perot(vars, R_fabry_perot_value)]


def calculate_beta_and_R(Reflection, Transmission):
    omega_beta_R_List = []
    for r, t in zip(Reflection, Transmission):
        beta, R = fsolve(equations, initial_guess, args=(t[1], r[1]))
        omega_beta_R_List.append([r[0], beta, R])
        if equations((beta, R), t[1], r[1])[0] + equations((beta, R), t[1], r[1])[1] > 0.001:
            print('Für beta = ' + str(beta) + ' und R = ' + str(R) + ' ist die Lsg: ' + str(equations((beta, R), t[1], r[1])))

    return omega_beta_R_List


reflection = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_GaAs_doped_res4_N50_normalized.DPT')
transmission = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\GaAs_doped_res4_N50_normalized.DPT')


beta = []
beta2 = []
for i in calculate_beta_and_R(bügeln(reflection, 1), bügeln(transmission, 1)):
    beta.append([i[0], i[1]])
plot_data(beta)
#for i in calculate_beta_and_R(reflection,transmission):
#    beta2.append([i[0], i[1]])
#plot_data(bügeln(beta2, 0.8))
save_and_open("foo")




def plotKomischeFunktion():
    read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_doped_res03_N50_new_normalized.DPT')
    
    komischeFunktion = []
    for s, r in zip(dataSample, dataReference):
        komischeFunktion.append([s[0], s[1] / r[1]])
    
    plot_data()
    plt.xlabel('wave number k / ' + r'$cm^{-1}$')    
    save_and_open("foo")
    
#plotKomischeFunktion()
    
    