from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file

filepathNormalized = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_normalized.DPT'
filepathSample = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_sample.DPT'
filepathReference = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_reference.DPT'



dataSample = read_dpt_file(filepathSample)
dataReference = read_dpt_file(filepathReference)



#zum Überprüfen, ob die normalized Daten gleich meiner berechnnung sind
#divided_data=[]
#for s, r in zip(dataSample, dataReference):
#    divided_data.append([s[0], s[1] / r[1]])
#plot_data(divided_data)

def plotGasAbs():
    plot_data(read_dpt_file(filepathNormalized))
    plt.xlabel('wave number \nu / ' + r'$cm^{-1}$', fontsize=17)
    plt.ylabel('Transmission T', fontsize=17)
    plt.tick_params(axis='both', which='major', labelsize=17)


    plt.text(2100, 0.62, r'$(1):\ \text{CO}_2$', fontsize=15, color='black')
    plt.text(3500, 0.68, r'$(2):\ \text{CO}_2, \text{H}_2\text{O}$', fontsize=15, color='black')
    plt.text(5100, 0.95, r'$(3):\ \text{H}_2\text{O}$', fontsize=15, color='black')
    plt.text(7000, 0.98, r'$(4):\ \text{H}_2\text{O}$', fontsize=15, color='black')

    save_and_open("GasAbsorption")
    
    
def plotSample():
    plot_data(read_dpt_file(filepathSample))
    plt.xlabel('wave number \nu / ' + r'$cm^{-1}$', fontsize=17)
    plt.ylabel('Transmission T / arbitrary units', fontsize=17)
    plt.tick_params(axis='both', which='major', labelsize=17)
    save_and_open("GasAbsorptionSample")



plotSample()
plotGasAbs()