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
from countTo50 import plot_Discrete_ns_for_comparison


e = const.e  # Elementary charge in Coulombs
hbar = const.hbar  # Reduced Planck's constant in J.s
h = const.h  # Planck's constant in J.s
m = const.m_e  # Electron mass in kg
epsilon_0 = const.epsilon_0  # Vacuum permittivity in F/m
c= const.c  # Speed of light in m/s

#moch willkürliche Werte
mu_GaAs = 0.036872 * const.m_e 
mu_GaAs_Lukas = 0.05607 * const.m_e
mu_GaSb = 0.0053207 * const.m_e
mu_GaSb_Lukas = 0.037188208 * const.m_e
m_e = const.m_e  # Electron mass in kg

factor_GaAs = e**4*8/(epsilon_0**2*m_e**4*c**3*h**5) * mu_GaAs_Lukas**3 * (2*np.pi*c)**4
factor_GaSb = e**4*8/(epsilon_0**2*m_e**4*c**3*h**5) * mu_GaSb_Lukas**3 * (2*np.pi*c)**4



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

samplesOhneSiUnVergleich = [
    (reflectionSiUnDo, transmissionSiUnDo, "Si Undoped", 530*10**-6,'blue'),
        (reflectionGaAsUnDo, transmissionGaAsUnDo, "GaAs Undoped", 470*10**-6, 'green'),
        (reflectionGaAsDo, transmissionGaAsDo, "GaAs Doped", 440*10**-6, 'red')
        
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
    nu_beta_R_List = []
    for r, t in zip(bügeln(ReflectionData,1), bügeln(TransmissionData,1)):
        initial_guess = nu_beta_R_List[-1][1:] if nu_beta_R_List else [10000, 0.3]  # Use last beta and R or default
        #initial_guess = [20000,0.3]
        beta, R = fsolve(equations, initial_guess, args=(t[1], r[1], d))
        if beta > 28000:
            beta = nu_beta_R_List[-1][1] if nu_beta_R_List else 10000  # Use last beta or default
        if R >1:
            R = nu_beta_R_List[-1][2] if nu_beta_R_List else 0.3 # R > 1 unphysical
        nu_beta_R_List.append([r[0], beta, R])
    return nu_beta_R_List


def calculate_kappa(ReflectionData, TransmissionData, d):
    nu_beta_R_List = calculate_nu_beta_and_R(ReflectionData, TransmissionData, d)
    nu_kappa_R_List = []
    for i in nu_beta_R_List:
        nu_kappa_R_List.append([i[0], i[1]/(2*2*np.pi*i[0]*100), i[2]]) # 100 wegen cm^-1
    return nu_kappa_R_List

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
    plt.ylabel('refractive index n')
    plt.legend(fontsize=12)
    save_and_open('RefractiveIndicesSmooth')
    
    
def plot_the_ns_vergleich(glättwert=0.02):
    plt.figure()
    for reflection, transmission, label, d, color in samplesOhneSiUnVergleich:
        n = []
        for i in calculate_k_n1_kappa(reflection, transmission, d):
            n.append([i[0], i[1]])
        plot_data(bügeln(n, glättwert), label=label, color = color)
    plt.xlabel(r'wave number $\nu$  / ' + r'$cm^{-1}$')
    plt.ylabel('refractive index n')
    plot_Discrete_ns_for_comparison()
    plt.ylim(2.8,4.2)
    plt.xlim(1800, 10500)
    
    plt.legend(fontsize=12)
    save_and_open('RefractiveIndicesVergleich')


     
def plot_the_betas(title="foo", glättwert=1):
    plt.figure()
    for reflection, transmission, label, d in samplesOhneSiUn:
        beta = []
        for i in calculate_nu_beta_and_R(bügeln(reflection, glättwert), bügeln(transmission, glättwert), d):
            beta.append([i[0], i[1]/100])
        plot_data(bügeln(beta, 0.05), label=label)
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
        
        kappa = [[i[0], i[1] * 10**3] for i in kappa]
        
        
        plot_data(bügeln(kappa, 0.05), label=label)

    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plt.ylabel(r'extinction coefficient $\kappa$ / ' + r'$ 10^{-3}\, cm^{-1}$')
    plt.legend(fontsize=12)
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
        komischeFunktion.append([nu, (epsilon_2*c**2*nu**2*4*np.pi**2*(100**2))**2])
    return komischeFunktion


def plotKomischeFunktion(reflection, transmission, k_min, k_max, k_min_regression, k_max_regression, d, factor, glättwert=0.01, title="foo", varFürFehler=15):
    komischeFunktion = bügeln(calculate_KomischeFunktion(reflection, transmission, d),glättwert)
    
    #Lin Reg Teil
    chopped = k_chopper(komischeFunktion, k_min=k_min_regression, k_max=k_max_regression)
    steig, x_intercept, y_intercept = linear_regression(chopped)
    x_vals = np.linspace(x_intercept, k_max_regression, 10)  # Generate x values for the line
    y_vals = steig * x_vals + y_intercept       # Calculate corresponding y values
    plt.plot(x_vals, y_vals*10**-56, label='Linear Fit', color='red', linestyle='-')  # Plot the line
    plt.legend(fontsize=18)
    
    # Fehler LinReg Teil
    steig1, x_intercept1, y_intercept1 = linear_regression(k_chopper(komischeFunktion,
                                                k_min=k_min_regression - varFürFehler, k_max=k_max_regression + varFürFehler))
    steig2, x_intercept2, y_intercept2 = linear_regression(k_chopper(komischeFunktion,
                                                k_min=k_min_regression - varFürFehler, k_max=k_max_regression - varFürFehler))
    steig3, x_intercept3, y_intercept3 = linear_regression(k_chopper(komischeFunktion,
                                                k_min=k_min_regression + varFürFehler, k_max=k_max_regression + varFürFehler))
    steig4, x_intercept4, y_intercept4 = linear_regression(k_chopper(komischeFunktion,
                                                k_min=k_min_regression + varFürFehler, k_max=k_max_regression - varFürFehler))
    
    steig_fehler = max(abs(steig1-steig), abs(steig2-steig), abs(steig3-steig), abs(steig4-steig))
    x_intercept_fehler = max(abs(x_intercept1-x_intercept), abs(x_intercept2-x_intercept), abs(x_intercept3-x_intercept), abs(x_intercept4-x_intercept))
    
    komischeFunktion56 = [[point[0], point[1] * 10**-56] for point in komischeFunktion]
    plot_data(k_chopper(komischeFunktion56,k_min, k_max), label=title)
    
    plt.axvline(x=k_min_regression, color='red', linestyle='--', label='Boundary for fit', linewidth=0.5)
    plt.legend(fontsize=16)
    plt.axvline(x=k_max_regression, color='red', linestyle='--', linewidth =0.5)
    
    
    steigKorr = steig/factor
    steigKorr_fehler = steig_fehler/factor
    pulseMatrixElement = steigKorr**(1/4)
    pulseMatrixElement_fehler = steigKorr_fehler / steigKorr /4 * pulseMatrixElement
    
    writeLatexMacro('bandgap_' + title.replace(' ','_'), x_intercept*100*c*hbar/e*2*np.pi, 'eV', x_intercept_fehler*100*c*hbar/e*2*np.pi)
    
    writeLatexMacro('pulseMatrixElement_' + title.replace(' ','_'), pulseMatrixElement, r'$m\, kg s^{-1}$', pulseMatrixElement_fehler)
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=19)   
    plt.ylabel(r'($\epsilon ^{\prime \prime} \omega^2)^2\, / \, 10^{56} \left[\frac{A}{Vms}\right]^2$', fontsize=19) 
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
        beta = i[1]+100 #make sure beta is posivive (constant shift of 1), negative beta is not physical
        komischeFunktion.append([nu, np.sqrt(beta)*nu*c*2*np.pi*100 * 10**-17])
    return komischeFunktion

def plotKomischeIndirectFunction(reflection, transmission, d, k_min, k_max, k1_min_regression, 
                                 k1_max_regression, k2_min_regression, k2_max_regression, glättwert=0.01, title="foo",
                                 varFürFehler=15):
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
    steigVerschoben1, x_interceptVerschoben1, y_interceptVerschoben1 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k1_min_regression - varFürFehler, k_max=k1_max_regression + varFürFehler))
    steigVerschoben2, x_interceptVerschoben2, y_interceptVerschoben2 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k1_min_regression - varFürFehler, k_max=k1_max_regression - varFürFehler))
    steigVerschoben3, x_interceptVerschoben3, y_interceptVerschoben3 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k1_min_regression + varFürFehler, k_max=k1_max_regression + varFürFehler))
    steigVerschoben4, x_interceptVerschoben4, y_interceptVerschoben4 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k1_min_regression + varFürFehler, k_max=k1_max_regression - varFürFehler))
    x_intercept_fehler1 = max(abs(x_interceptVerschoben1-x_intercept1), abs(x_interceptVerschoben2-x_intercept1),
                             abs(x_interceptVerschoben3-x_intercept1), abs(x_interceptVerschoben4-x_intercept1))
    
    
    
    #Lin Reg Teil 2
    chopped = k_chopper(komischeIndirFunktion, k_min=k2_min_regression, k_max=k2_max_regression)
    steig2, x_intercept2, y_intercept2 = linear_regression(chopped)
    x_vals = np.linspace(x_intercept2, k2_max_regression+100, 10)  # Generate x values for the line
    y_vals = steig2 * x_vals + y_intercept2       # Calculate corresponding y values
    plt.plot(x_vals, y_vals, label='Linear Fit', color='red', linestyle='-')  # Plot the line
    plt.axvline(x=k2_min_regression, color='red', linestyle='--', linewidth=1)
    plt.axvline(x=k2_max_regression, color='red', linestyle='--', linewidth =1)
    plt.legend()
    steig2Verschoben1, x_intercept2Verschoben1, y_intercept2Verschoben1 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k2_min_regression - varFürFehler, k_max=k2_max_regression + varFürFehler))
    steig2Verschoben2, x_intercept2Verschoben2, y_intercept2Verschoben2 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k2_min_regression - varFürFehler, k_max=k2_max_regression - varFürFehler))
    steig2Verschoben3, x_intercept2Verschoben3, y_intercept2Verschoben3 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k2_min_regression + varFürFehler, k_max=k2_max_regression + varFürFehler))
    steig2Verschoben4, x_intercept2Verschoben4, y_intercept2Verschoben4 = linear_regression(k_chopper(komischeIndirFunktion,
                                                k_min=k2_min_regression + varFürFehler, k_max=k2_max_regression - varFürFehler))
    x_intercept_fehler2 = max(abs(x_interceptVerschoben1-x_intercept1), abs(x_interceptVerschoben2-x_intercept1),
                             abs(x_interceptVerschoben3-x_intercept1), abs(x_interceptVerschoben4-x_intercept1))

    x_fehler_gesamt = (x_intercept_fehler1 + x_intercept_fehler2) /2
    
    bandgap = (x_intercept1 + x_intercept2)/2*100*c*hbar/e*2*np.pi
    bandgap_fehler = x_fehler_gesamt *100*c*hbar/e*2*np.pi
    
    hquerOMEGA = (x_intercept2-x_intercept1)/2*100*c*hbar/e*2*np.pi
    hquerOMEGA_fehler = x_fehler_gesamt*100*c*hbar/e*2*np.pi
    
    writeLatexMacro("OMEGA_" + title, hquerOMEGA/hbar, 'Hz', hquerOMEGA_fehler/hbar) # ist der gleiche fehler wie für die bandgap
     
    writeLatexMacro("bandgap_" + title, bandgap, 'eV', bandgap_fehler)
    writeLatexMacro("hquerOMEGA_" + title, hquerOMEGA, 'eV', hquerOMEGA_fehler) # ist der gleiche fehler wie für die bandgap
    
    plt.legend(fontsize=18)
    
    plot_data(k_chopper(komischeIndirFunktion, k_min=k_min, k_max=k_max), label=title)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=18)
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    plt.ylabel(r'$\omega \cdot \sqrt{\beta}$ / ' + r'$10^{17}\,s^{-1}\,m^{-1/2}$', fontsize=18)
    save_and_open(filename=title)

def plot_reflections(glättwert=0.01):
    plt.figure()
    
    for reflection, transmission, label, d in samplesOhneSiUn:
        plot_data(bügeln(reflection, glättwert), label=label)
    plt.ylabel('Reflection R')
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
    plot_data(bügeln(reflectionSiDo, glättwert), label = "Si Doped") 
    
    plt.legend(loc='upper right', bbox_to_anchor=(0.53, 0.73))
    #plt.legend()
    save_and_open("Low_Res_Reflections")
    
def plot_transmissions(glättwert=0.01):
    plt.figure()
    
    for reflection, transmission, label, d in samplesOhneSiUn:
        plot_data(bügeln(transmission, glättwert), label=label)
        plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$')
        plt.ylabel('Transmission T')
    plot_data(bügeln(transmissionSiDo, glättwert), label="Si Doped")    
    plt.legend()
    save_and_open("Low_Res_Transmissions")






transmissionGaAsDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_doped_res03_N50_new_normalized.DPT')
#2000-9000
transmissionGaAsUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_undoped_res03_N50_new_normalized.DPT')
#2000-10500
transmissionGaSbDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaSb_doped_res03_N50_new_normalized.DPT')
#GaSb nicht möglich wegen rauschen
transmissionSiUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\Si_undoped_res03_N50_normalized.DPT')
#SiUnDo 2000-8500

samplesOhneSiUnDiscrete = [(transmissionSiUnDo, "Si Undoped", 530*10**-6, 0.0035, 2000, 8500),
        (transmissionGaAsUnDo, "GaAs Undoped", 470*10**-6, 0.004, 2000, 10500),
        (transmissionGaAsDo, "GaAs Doped", 440*10**-6, 0.004, 2000, 9000)
    ]

from countTo50 import calculate_n
def plot_epxilon_stich_with_scuffed_kappa(windowsize=200):
    colors = ['blue', 'green', 'red', 'purple', 'orange']  # Add a list of colors
    for idx, (data, label, d, prominence, nu_min, nu_max) in enumerate(samplesOhneSiUnDiscrete):
        nus = np.array(range(nu_min, nu_max, 500))
        ns = []
        nerrors =[]
        kappas = []
        for nu in nus:
            ns.append(calculate_n(data, d, nu, windowsize, prominence=prominence))
            nerrors.append(0.02 * calculate_n(data, d, nu, windowsize, prominence=prominence))
        
        for reflection, transmission, label_kontinuierlich, d,_ in samplesOhneSiUnVergleich:
            if label == label_kontinuierlich:
                kappa_liste = np.array(calculate_kappa(reflection, transmission, d))
                                
                kappa_liste = sorted(kappa_liste, key=lambda x: x[0])
                kappas = np.interp(nus, [item[0] for item in kappa_liste], [item[1] for item in kappa_liste])
        
        
        epsilon_striche = [n**2-kappa**2 for (n,kappa) in zip(ns, kappas)]
        epsilon_strich_fehler = [0.02 * np.sqrt(2)/2 * epsilon_strich for epsilon_strich in epsilon_striche]
        #epsilon_2_strich = [2*n*kappa for (n,kappa) in zip(ns, kappas)]
        #epsilon_2_strich_fehler = [0.02 * epsilon_2_strich for epsilon_2_strich in epsilon_2_striche]

        plt.scatter(nus, epsilon_striche, color=colors[idx], label=label)
        plt.errorbar(nus, epsilon_striche, yerr=epsilon_strich_fehler, fmt='o', color=colors[idx % len(colors)], capsize=5)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=18)
    plt.ylabel(r'$\epsilon^\prime$', fontsize=19)
    plt.grid()
    plt.legend(fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    import DataPlotter as plotter
    plotter.save_and_open(filename='epsilon1StrichDiskret')                
    

def plot_epxilon_2_stich_with_scuffed_kappa(windowsize=200):
    colors = ['blue', 'green', 'red', 'purple', 'orange']  # Add a list of colors
    for idx, (data, label, d, prominence, nu_min, nu_max) in enumerate(samplesOhneSiUnDiscrete):
        nus = np.array(range(nu_min, nu_max, 500))
        ns = []
        nerrors =[]
        kappas = []
        for nu in nus:
            ns.append(calculate_n(data, d, nu, windowsize, prominence=prominence))
            nerrors.append(0.02 * calculate_n(data, d, nu, windowsize, prominence=prominence))
        
        for reflection, transmission, label_kontinuierlich, d,_ in samplesOhneSiUnVergleich:
            if label == label_kontinuierlich:
                kappa_liste = np.array(calculate_kappa(reflection, transmission, d))
                                
                kappa_liste = sorted(kappa_liste, key=lambda x: x[0])
                kappas = np.interp(nus, [item[0] for item in kappa_liste], [item[1] for item in kappa_liste])
        
        
        #epsilon_striche = [n**2-kappa**2 for (n,kappa) in zip(ns, kappas)]
        #epsilon_strich_fehler = [0.02 * np.sqrt(2)/2 * epsilon_strich for epsilon_strich in epsilon_striche]
        epsilon_2_striche = [2*n*kappa for (n,kappa) in zip(ns, kappas)]
        epsilon_2_strich_fehler = [np.abs(0.02 * epsilon_2_strich) for epsilon_2_strich in epsilon_2_striche]

        plt.scatter(nus, epsilon_2_striche, color=colors[idx], label=label)
        plt.errorbar(nus, epsilon_2_striche, yerr=epsilon_2_strich_fehler, fmt='o', color=colors[idx % len(colors)], capsize=5)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=18)
    plt.ylabel(r'$\epsilon^{\prime\prime}$', fontsize=19)
    plt.grid()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=18)
    import DataPlotter as plotter
    plotter.save_and_open(filename='epsilon2StrichDiskret')   
                


#plotKomischeIndirectFunction(reflectionSiUnDo, transmissionSiUnDo, d=530*10**-6, k_min=7500, k_max=11000,
#                              k1_min_regression=8450, k1_max_regression=9100, k2_min_regression=9450, k2_max_regression=10200,
#                              glättwert=0.01, title="Si_Undoped")


#plot_the_kappas(glättwert=0.9, title="kappa")

#plot_the_betas(title="betas", glättwert=0.9)

#plot_the_ns(0.005)
plot_epxilon_2_stich_with_scuffed_kappa()
plot_epxilon_stich_with_scuffed_kappa()

#plot_the_ns_vergleich()

#Die macht keinen Sinn, ist indirekt!! plotKomischeFunktion(reflectionSiUnDo, transmissionSiUnDo, 9000, 11000, 10050, 10250, samplesOhneSiUn[0][3], glättwert=0.1, title="Si Undoped")
#plotKomischeFunktion(reflectionGaSbDo, transmissionGaSbDo, 5000, 6000, 5600, 5680, samplesOhneSiUn[1][3], factor_GaSb ,glättwert=0.1, title="GaSb Doped")
#plotKomischeFunktion(reflectionGaAsUnDo, transmissionGaAsUnDo, 11000, 11400, 11230, 11280, samplesOhneSiUn[2][3], factor_GaAs, glättwert=0.1, title="GaAs Undoped")
#plotKomischeFunktion(reflectionGaAsDo, transmissionGaAsDo, 10500, 12000, 11120, 11200, samplesOhneSiUn[3][3],factor_GaAs, glättwert=0.03, title="GaAs Doped")


#plot_reflections(glättwert=0.02)
#plot_transmissions(glättwert=1)




