from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
import numpy as np
import scipy.constants as const
from scipy.optimize import fsolve
from Frauen import bügeln


e = const.e  # Elementary charge in Coulombs
hbar = const.hbar  # Reduced Planck's constant in J.s
m = const.m_e  # Electron mass in kg
epsilon_0 = const.epsilon_0  # Vacuum permittivity in F/m
c= const.c  # Speed of light in m/s

#moch willkürliche Werte
mu_GaAs = 0.036872 * const.m_e  # Effective mass of electron in GaAs (müssen quelle finden!!)
m_e = const.m_e  # Electron mass in kg
factor = e*e*(2*mu_GaAs)**(3/2)*2*np.pi/(epsilon_0*m_e**2*hbar**3)    
d = 440*10**-6  # Thickness of the sample in m



#calculate n

initial_guess = [1,1] # Initial guess for beta and R

def T_fabry_perot(vars, T_fabry_perot_value):
    beta, R_halbraum = vars
    return T_fabry_perot_value - (1-R_halbraum)**2*np.exp(-beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d))

def R_fabry_perot(vars, R_fabry_perot_value):
    beta, R_halbraum = vars
    return R_fabry_perot_value -(R_halbraum + (1-R_halbraum)**2*R_halbraum*np.exp(-2*beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d)))

def equations(vars, T_fabry_perot_value, R_fabry_perot_value):
    return [T_fabry_perot(vars, T_fabry_perot_value), R_fabry_perot(vars, R_fabry_perot_value)]


def calculate_beta_and_R(ReflectionData, TransmissionData):
    k_beta_R_List = []
    for r, t in zip(bügeln(ReflectionData,1), bügeln(TransmissionData,1)):
        initial_guess = k_beta_R_List[-1][1:] if k_beta_R_List else [10000, 0.3]  # Use last beta and R or default
        #initial_guess = [20000,0.3]
        beta, R = fsolve(equations, initial_guess, args=(t[1], r[1]))
        if beta > 28000:
            beta = k_beta_R_List[-1][1] if k_beta_R_List else 10000  # Use last beta or default
        if R >1:
            R = k_beta_R_List[-1][2] if k_beta_R_List else 0.3 
        k_beta_R_List.append([r[0], beta, R])
    return k_beta_R_List


def calculate_kappa(ReflectionData, TransmissionData):
    k_beta_R_List = calculate_beta_and_R(ReflectionData, TransmissionData)
    k_kappa_R_List = []
    for i in k_beta_R_List:
        k_kappa_R_List.append([i[0], i[1]/(2*i[0]*100), i[2]]) # 100 wegen cm^-1
    return k_kappa_R_List

def calculate_k_n1_kappa(ReflectionData, TransmissionData):
    k_kappa_R_List = calculate_kappa(ReflectionData, TransmissionData)
    n = []
    for i in k_kappa_R_List:
        k = i[1]
        R = i[2]
        n1 = (1+R+np.sqrt(-k**2-R**2*k**2+2*R*(2+k**2)))/(1-R)
        n2 = (-1-R+np.sqrt(-k**2-R**2*k**2+2*R*(2+k**2)))/(R-1)
                
        n.append([i[0], n1,i[1]])
    return n

def plot_the_ns(glättwert=0.02):
    plt.figure()
    for reflection, transmission, label in samplesOhneSiUn:
        n = []
        for i in calculate_k_n1_kappa(reflection, transmission):
            n.append([i[0], i[1]])
        plot_data(bügeln(n, glättwert), label=label)
    plt.xlabel('wave number k / ' + r'$cm^{-1}$')
    plt.ylabel('n')
    save_and_open('RefractiveIndicesSmooth')
    
    
def plot_the_betas():
    plt.figure()
    for reflection, transmission, label in samplesOhneSiUn:
        beta = []
        for i in calculate_beta_and_R(reflection, transmission):
            beta.append([i[0], i[1]])
        plot_data(beta, label=label)
        plt.xlabel('wave number k / ' + r'$cm^{-1}$')
        plt.ylabel('beta')

    plt.xlabel('wave number k / ' + r'$cm^{-1}$')
    plt.ylabel('beta')
    save_and_open('foo')
    
    
def calculate_k_epsilon_2(reflection, transmission):
    k_n_kappa = calculate_k_n1_kappa(reflection, transmission)
    epsilon_2 = []
    for i in k_n_kappa:
        n = i[1]
        kappa = i[2]
        epsilon_2.append([i[0], 2*n*kappa])
    return epsilon_2


def k_chopper(data, k_min, k_max):
    return [point for point in data if k_min <= point[0] <= k_max]


def calculate_KomischeFunktion(reflection, transmission, glättwert):
    
    k_epsilon_2 = calculate_k_epsilon_2(reflection, transmission)
    
    komischeFunktion = []
    for k, epsilon_2 in k_epsilon_2:
        komischeFunktion.append([k, (epsilon_2*c**2*k**2)**2])
    return komischeFunktion


def plotKomischeFunktion(reflection, transmission, glättwert=0.01):
    komischeFunktion = calculate_KomischeFunktion(reflection, transmission, glättwert)
    plot_data(k_chopper(bügeln(komischeFunktion, glättwert),9000,13000))
    
    chopped = k_chopper(calculate_KomischeFunktion(reflectionGaAsDo, transmissionGaAsDo, glättwert=0.01), 11000, 11500)
    linear_regression(chopped)
    
    
    
    
    
    plt.xlabel('wave number k / ' + r'$cm^{-1}$')    
    save_and_open("foo")
    
def linear_regression(komischeFunktion):
    """
    Perform a linear regression on the given x and y data.

    Parameters:
        x (list or np.ndarray): Independent variable data.
        y (list or np.ndarray): Dependent variable data.

    Returns:
        tuple: (slope, intercept, x_intercept) of the best-fit line.
    """
    x = np.array([point[0] for point in komischeFunktion])
    y = np.array([point[1] for point in komischeFunktion])
    x = np.array(komischeFunktion[0])
    y = np.array(komischeFunktion[1])
    n = len(x)
    slope = (n * np.sum(x * y) - np.sum(x) * np.sum(y)) / (n * np.sum(x**2) - np.sum(x)**2)
    intercept = (np.sum(y) - slope * np.sum(x)) / n
    x_intercept = -intercept / slope if slope != 0 else None
    return slope, intercept, x_intercept

    


reflectionGaAsDo = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_GaAs_doped_res4_N50_normalized.DPT')
transmissionGaAsDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\GaAs_doped_res4_N50_normalized.DPT')

reflectionGaAsUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_GaAs_undoped_res4_N50_normalized.DPT')
transmissionGaAsUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\GaAs_undoped_res4_N50_normalized.DPT')

reflectionGaSbDo = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_GaSb_doped_res4_N50_normalized.DPT')
transmissionGaSbDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\GaSb_doped_res4_N50_normalized.DPT')

reflectionSiUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_Si_undoped_res4_N50_normalized.DPT')
transmissionSiUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\Si_undoped_res4_N50_normalized.DPT')

reflectionSiDo = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_Si_doped_res4_N50_normalized.DPT')
transmissionSiDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\Si_doped_res4_N50_normalized.DPT')

samplesOhneSiUn = [
        (reflectionGaAsDo, transmissionGaAsDo, "GaAs Doped"),
        (reflectionGaAsUnDo, transmissionGaAsUnDo, "GaAs Undoped"),
        (reflectionGaSbDo, transmissionGaSbDo, "GaSb Doped"),
        (reflectionSiUnDo, transmissionSiUnDo, "Si UnDoped")
    ]

# n1=[]
# for i in calculate_n(reflection, transmission):
#     n1.append([i[0], i[1]])
# plot_data(n1)
# save_and_open("foo")



#plot_the_betas()

# kappa =[]
# for i in calculate_kappa(reflectionGaAsDo, transmissionGaAsDo):
#     kappa.append([i[0], i[1]])
# plot_data(kappa)
# plt.xlabel('wave number k / ' + r'$cm^{-1}$')
# plt.ylabel('kappa')
# save_and_open('foo')


#plotKomischeFunktion(reflectionGaAsDo, transmissionGaAsDo, glättwert=0.01)
#plot_the_ns(0.005)




