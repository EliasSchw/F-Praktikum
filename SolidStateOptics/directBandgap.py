from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
import numpy as np
import scipy.constants as const

e = const.e  # Elementary charge in Coulombs
hbar = const.hbar  # Reduced Planck's constant in J.s
m = const.m_e  # Electron mass in kg
epsilon_0 = const.epsilon_0  # Vacuum permittivity in F/m

#moch willkürliche Werte
mu_GaAs = 0.036872 * const.m_e  # Effective mass of electron in GaAs (müssen quelle finden!!)
m_e = const.m_e  # Electron mass in kg
factor = e*e*(2*mu_GaAs)**(3/2)*2*np.pi/(epsilon_0*m_e**2*hbar**3)    





def plotKomischeFunktion():
    read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_doped_res03_N50_new_normalized.DPT')
    
    komischeFunktion = []
    for s, r in zip(dataSample, dataReference):
        komischeFunktion.append([s[0], s[1] / r[1]])
    
    plot_data()
    plt.xlabel('wave number k / ' + r'$cm^{-1}$')    
    save_and_open("foo")
    
plotKomischeFunktion()
    
    