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


plot_data(read_dpt_file(filepathNormalized))
plt.xlabel('wave number k / ' + r'$cm^{-1}$')
plt.ylabel('Transmission T')

save_and_open("GasAbsorption")