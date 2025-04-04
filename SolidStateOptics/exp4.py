from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
import numpy as np
import scipy.constants as const
from scipy.optimize import fsolve
from Frauen import bügeln
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from macroswriter import writeLatexMacro


e = const.e  # Elementary charge in Coulombs
hbar = const.hbar  # Reduced Planck's constant in J.s
m = const.m_e  # Electron mass in kg
epsilon_0 = const.epsilon_0  # Vacuum permittivity in F/m
c= const.c  # Speed of light in m/s

#moch willkürliche Werte
mu_GaAs = 0.036872 * const.m_e  # Effective mass of electron in GaAs (müssen quelle finden!!)
m_e = const.m_e  # Electron mass in kg
factor = e*e*(2*mu_GaAs)**(3/2)*2*np.pi/(epsilon_0*m_e**2*hbar**3)    

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
        (reflectionGaAsDo, transmissionGaAsDo, "GaAs Doped", 440*10**-6),
        (reflectionGaAsUnDo, transmissionGaAsUnDo, "GaAs Undoped", 470*10**-6),
        (reflectionGaSbDo, transmissionGaSbDo, "GaSb Doped", 500*10**-6),
        (reflectionSiUnDo, transmissionSiUnDo, "Si Undoped", 530*10**-6)
    ]

def T_fabry_perot(vars, T_fabry_perot_value, d):
    beta, R_halbraum = vars
    return T_fabry_perot_value - (1-R_halbraum)**2*np.exp(-beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d))

def R_fabry_perot(vars, R_fabry_perot_value, d):
    beta, R_halbraum = vars
    return R_fabry_perot_value -(R_halbraum + (1-R_halbraum)**2*R_halbraum*np.exp(-2*beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d)))

def equations(vars, T_fabry_perot_value, R_fabry_perot_value, d):
    return [T_fabry_perot(vars, T_fabry_perot_value, d), R_fabry_perot(vars, R_fabry_perot_value, d)]


def calculate_beta_and_R(ReflectionData, TransmissionData, d):
    k_beta_R_List = []
    for r, t in zip(bügeln(ReflectionData,1), bügeln(TransmissionData,1)):
        initial_guess = k_beta_R_List[-1][1:] if k_beta_R_List else [10000, 0.3]  # Use last beta and R or default
        #initial_guess = [20000,0.3]
        beta, R = fsolve(equations, initial_guess, args=(t[1], r[1], d))
        if beta > 28000:
            beta = k_beta_R_List[-1][1] if k_beta_R_List else 10000  # Use last beta or default
        if R >1:
            R = k_beta_R_List[-1][2] if k_beta_R_List else 0.3 
        k_beta_R_List.append([r[0], beta, R])
    return k_beta_R_List


def calculate_kappa(ReflectionData, TransmissionData, d):
    k_beta_R_List = calculate_beta_and_R(ReflectionData, TransmissionData, d)
    k_kappa_R_List = []
    for i in k_beta_R_List:
        k_kappa_R_List.append([i[0], i[1]/(2*i[0]*100), i[2]]) # 100 wegen cm^-1
    return k_kappa_R_List

def calculate_k_n1_kappa(ReflectionData, TransmissionData, d):
    k_kappa_R_List = calculate_kappa(ReflectionData, TransmissionData, d)
    n = []
    for i in k_kappa_R_List:
        k = i[1]
        R = i[2]
        n1 = (1+R+np.sqrt(-k**2-R**2*k**2+2*R*(2+k**2)))/(1-R)
        #n2 = (-1-R+np.sqrt(-k**2-R**2*k**2+2*R*(2+k**2)))/(R-1)
                
        n.append([i[0], n1,i[1]])
    return n

def plot_the_ns(glättwert=0.02):
    plt.figure()
    for reflection, transmission, label, d in samplesOhneSiUn:
        n = []
        for i in calculate_k_n1_kappa(reflection, transmission, d):
            n.append([i[0], i[1]])
        plot_data(bügeln(n, glättwert), label=label)
    plt.xlabel(r'wave number $\nu$  / ' + r'$cm^{-1}$')
    plt.ylabel('n')
    save_and_open('RefractiveIndicesSmooth')
    
    
def plot_the_betas(title="foo"):
    plt.figure()
    for reflection, transmission, label, d in samplesOhneSiUn:
        beta = []
        for i in calculate_beta_and_R(reflection, transmission, d):
            beta.append([i[0], i[1]/100])
        plot_data(beta, label=label)
        plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
        plt.ylabel('beta')

    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plt.ylabel(r'absorption coefficient $\beta$ / ' + r'$cm^{-1}$')
    plt.legend()
    save_and_open(filename=title)
    
    
def calculate_k_epsilon_2(reflection, transmission, d):
    k_n_kappa = calculate_k_n1_kappa(reflection, transmission, d)
    epsilon_2 = []
    for i in k_n_kappa:
        n = i[1]
        kappa = i[2]
        epsilon_2.append([i[0], 2*n*kappa])
    return epsilon_2


def k_chopper(data, k_min, k_max):
    return [point for point in data if k_min <= point[0] <= k_max]


def calculate_KomischeFunktion(reflection, transmission, d):
    
    k_epsilon_2 = calculate_k_epsilon_2(reflection, transmission, d)
    
    komischeFunktion = []
    for k, epsilon_2 in k_epsilon_2:
        komischeFunktion.append([k, (epsilon_2*c**2*k**2)**2])
    return komischeFunktion


def plotKomischeFunktion(reflection, transmission, k_min, k_max, k_min_regression, k_max_regression, d, glättwert=0.01, title="foo"):
    komischeFunktion = bügeln(calculate_KomischeFunktion(reflection, transmission, d),glättwert)
    
    chopped = k_chopper(komischeFunktion, k_min=k_min_regression, k_max=k_max_regression)
    steig, x_intercept, y_intercept = linear_regression(chopped)
    
    x_vals = np.linspace(x_intercept, k_max_regression, 10)  # Generate x values for the line
    y_vals = steig * x_vals + y_intercept       # Calculate corresponding y values
    plt.plot(x_vals, y_vals, label='Linear Fit', color='red', linestyle='-')  # Plot the line
    plt.legend()
        
    plot_data(k_chopper(komischeFunktion,k_min, k_max), label=title)
    
    # Plot a vertical line at k_min_regression
    plt.axvline(x=k_min_regression, color='blue', linestyle='--', label='boundary for fit', linewidth=0.5)
    plt.legend()
    plt.axvline(x=k_max_regression, color='blue', linestyle='--', linewidth =0.5)
    
    writeLatexMacro('bandgap_' + title, x_intercept*100*c*hbar/e*2*np.pi, 'eV')
    
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')    
    save_and_open(filename=title)
    
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
    # Removed incorrect overwriting of x and y
    n = len(x)
    slope = (n * np.sum(x * y) - np.sum(x) * np.sum(y)) / (n * np.sum(x**2) - np.sum(x)**2)
    y_intercept = (np.sum(y) - slope * np.sum(x)) / n
    x_intercept = -y_intercept / slope if slope != 0 else None
    return slope, x_intercept, y_intercept

    





# kappa =[]
# for i in calculate_kappa(reflectionGaAsDo, transmissionGaAsDo):
#     kappa.append([i[0], i[1]])
# plot_data(kappa)
# plt.xlabel('wave number k / ' + r'$cm^{-1}$')
# plt.ylabel('kappa')
# save_and_open('foo')












plotKomischeFunktion(reflectionGaAsDo, transmissionGaAsDo, 10500, 12000, 11120, 11200, samplesOhneSiUn[0][3], glättwert=0.03, title="GaAs Doped")
plotKomischeFunktion(reflectionGaAsUnDo, transmissionGaAsUnDo, 11000, 11400, 11230, 11280, samplesOhneSiUn[1][3], glättwert=0.1, title="GaAs Undoped")
plotKomischeFunktion(reflectionGaSbDo, transmissionGaSbDo, 5000, 6000, 5600, 5680, samplesOhneSiUn[2][3], glättwert=0.1, title="GaSb Doped")
plotKomischeFunktion(reflectionSiUnDo, transmissionSiUnDo, 9000, 11000, 10050, 10250, samplesOhneSiUn[3][3], glättwert=0.1, title="Si Undoped")


plot_the_ns(0.005)

plot_the_betas()





