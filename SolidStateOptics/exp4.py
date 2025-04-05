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

d_SiDo = 500*10**-6

samplesOhneSiUn = [
    (reflectionSiUnDo, transmissionSiUnDo, "Si Undoped", 530*10**-6),
        (reflectionGaSbDo, transmissionGaSbDo, "GaSb Doped", 500*10**-6),
        (reflectionGaAsUnDo, transmissionGaAsUnDo, "GaAs Undoped", 470*10**-6),
        (reflectionGaAsDo, transmissionGaAsDo, "GaAs Doped", 440*10**-6)
        
    ]

def T_fabry_perot(vars, T_fabry_perot_value, d):
    beta, R_halbraum = vars
    return T_fabry_perot_value - (1-R_halbraum)**2*np.exp(-beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d))

def R_fabry_perot(vars, R_fabry_perot_value, d):
    beta, R_halbraum = vars
    return R_fabry_perot_value -(R_halbraum + (1-R_halbraum)**2*R_halbraum*np.exp(-2*beta*d)/(1-R_halbraum**2*np.exp(-2*beta*d)))

def equations(vars, T_fabry_perot_value, R_fabry_perot_value, d):
    return [T_fabry_perot(vars, T_fabry_perot_value, d), R_fabry_perot(vars, R_fabry_perot_value, d)]


def calculate_nu_beta_and_R(ReflectionData, TransmissionData, d):
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
    k_beta_R_List = calculate_nu_beta_and_R(ReflectionData, TransmissionData, d)
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
    
    
def plot_the_betas(title="foo", glättwert=1):
    plt.figure()
    for reflection, transmission, label, d in samplesOhneSiUn:
        beta = []
        for i in calculate_nu_beta_and_R(bügeln(reflection, glättwert), bügeln(transmission, glättwert), d):
            beta.append([i[0], i[1]/100])
        plot_data(beta, label=label)
        plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
        plt.ylabel('beta')

    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plt.ylabel(r'absorption coefficient $\beta$ / ' + r'$cm^{-1}$')
    plt.legend()
    save_and_open(filename=title)
    
    
def plot_the_kappas(title="foo", glättwert = 1):
    plt.figure()
    for reflection, transmission, label, d in samplesOhneSiUn:
        kappa = []
        for i in calculate_kappa(bügeln(reflection, glättwert), bügeln(transmission,glättwert), d):
            kappa.append([i[0], i[1]])
        plot_data(kappa, label=label)
        plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
        plt.ylabel('kappa')

    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plt.ylabel(r'absorption coefficient $\kappa$ / ' + r'$cm^{-1}$')
    plt.legend()
    save_and_open(filename=title)
    
    
def calculate_nu_epsilon_2(reflection, transmission, d):
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
    
    nu_epsilon_2 = calculate_nu_epsilon_2(reflection, transmission, d)
    
    komischeFunktion = []
    for nu, epsilon_2 in nu_epsilon_2:
        komischeFunktion.append([nu, (epsilon_2*c**2*nu**2*4*np.pi**2/(100**2))**2])
    return komischeFunktion


def plotKomischeFunktion(reflection, transmission, k_min, k_max, k_min_regression, k_max_regression, d, glättwert=0.01, title="foo"):
    komischeFunktion = bügeln(calculate_KomischeFunktion(reflection, transmission, d),glättwert)
    
    #Lin Reg Teil
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
    
    pulseMatrixElement = np.sqrt(steig/factor)
    writeLatexMacro('pulseMatrixElement_' + title.replace(' ','_'), pulseMatrixElement, '??')
    
    writeLatexMacro('bandgap_' + title.replace(' ','_'), x_intercept*100*c*hbar/e*2*np.pi, 'eV')
    
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

    
def calculate_komische_indirect_function(reflection, transmission, d):    
    komischeFunktion = []
    for i in calculate_nu_beta_and_R(reflection, transmission, d):
        nu = i[0]
        beta = i[1]/100+1 #convert to cm^-1 and make sure beta is posivive (constant shift of 1), negative beta is not physical
        komischeFunktion.append([nu, np.sqrt(beta)*nu*c*2*np.pi])
    return komischeFunktion

def plotKomischeIndirectFunction(reflection, transmission, d, k_min, k_max, k1_min_regression, 
                                 k1_max_regression, k2_min_regression, k2_max_regression, glättwert=0.01, title="foo"):
    komischeIndirFunktion = bügeln(calculate_komische_indirect_function(reflection, transmission, d), glättwert)
    
    
    #Lin Reg Teil 1
    chopped = k_chopper(komischeIndirFunktion, k_min=k1_min_regression, k_max=k1_max_regression)
    steig1, x_intercept1, y_intercept1 = linear_regression(chopped)
    x_vals = np.linspace(x_intercept1, k1_max_regression+100, 10)  # Generate x values for the line
    y_vals = steig1 * x_vals + y_intercept1       # Calculate corresponding y values
    plt.plot(x_vals, y_vals, label='Linear Fit', color='purple', linestyle='-')  # Plot the line
    plt.axvline(x=k1_min_regression, color='purple', linestyle='--', label='', linewidth=1)
    plt.axvline(x=k1_max_regression, color='purple', linestyle='--', linewidth =1)
    plt.legend()
    
    
    #Lin Reg Teil 2
    chopped = k_chopper(komischeIndirFunktion, k_min=k2_min_regression, k_max=k2_max_regression)
    steig2, x_intercept2, y_intercept2 = linear_regression(chopped)
    x_vals = np.linspace(x_intercept2, k2_max_regression+100, 10)  # Generate x values for the line
    y_vals = steig2 * x_vals + y_intercept2       # Calculate corresponding y values
    plt.plot(x_vals, y_vals, label='Linear Fit', color='red', linestyle='-')  # Plot the line
    plt.axvline(x=k2_min_regression, color='red', linestyle='--', linewidth=1)
    plt.axvline(x=k2_max_regression, color='red', linestyle='--', linewidth =1)
    plt.legend()
    
    bandgap = (x_intercept1+x_intercept2)/2*100*c*hbar/e*2*np.pi
    writeLatexMacro("bandgap_" + title, bandgap, 'eV')
    hquerOMEGA = (x_intercept2-x_intercept1)/2*100*c*hbar/e*2*np.pi
    writeLatexMacro("hquerOMEGA_" + title, hquerOMEGA, 'eV')
    
    plt.legend()
    
    plot_data(k_chopper(komischeIndirFunktion, k_min=k_min, k_max=k_max), label=title)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plt.ylabel(r'$ck \cdot \sqrt{\beta}$ / ' + r'??')
    save_and_open(filename=title)






# beta_data = calculate_k_beta_and_R(reflectionSiUnDo, transmissionSiUnDo, samplesOhneSiUn[0][3])
# beta_plot = [[point[0], point[1] / 100 + 1] for point in beta_data]  # Convert beta to cm^-1
# plt.figure()
# plot_data(beta_plot, label="Si Undoped")
# plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
# plt.ylabel(r'absorption coefficient $\beta$ / ' + r'$cm^{-1}$')
# plt.legend()
# save_and_open(filename="SiUnDo_beta_plot")


#plotKomischeIndirectFunction(reflectionSiUnDo, transmissionSiUnDo, d=530*10**-6, k_min=7500, k_max=11000,
#                              k1_min_regression=8450, k1_max_regression=9100, k2_min_regression=9450, k2_max_regression=10200,
#                              glättwert=0.01, title="Si_Undoped")


# plot_the_kappas(glättwert=0.9, title="kappa")

# plot_the_betas(title="betas", glättwert=0.9)

# plot_the_ns(0.005)

#plotKomischeFunktion(reflectionSiUnDo, transmissionSiUnDo, 9000, 11000, 10050, 10250, samplesOhneSiUn[0][3], glättwert=0.1, title="Si Undoped")
#plotKomischeFunktion(reflectionGaSbDo, transmissionGaSbDo, 5000, 6000, 5600, 5680, samplesOhneSiUn[1][3], glättwert=0.1, title="GaSb Doped")
#plotKomischeFunktion(reflectionGaAsUnDo, transmissionGaAsUnDo, 11000, 11400, 11230, 11280, samplesOhneSiUn[2][3], glättwert=0.1, title="GaAs Undoped")
#plotKomischeFunktion(reflectionGaAsDo, transmissionGaAsDo, 10500, 12000, 11120, 11200, samplesOhneSiUn[3][3], glättwert=0.03, title="GaAs Doped")






