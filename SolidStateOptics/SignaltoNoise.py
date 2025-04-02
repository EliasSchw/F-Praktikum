from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file

filepathNormalizedN10 = '/Users/lukashein/Desktop/F-Praktikum-2/SolidStateOptics/RawData/StoN/StoN_res4_N10_normalized.DPT'
filepathNormalizedN20 = '/Users/lukashein/Desktop/F-Praktikum-2/SolidStateOptics/RawData/StoN/StoN_res4_N20_normalized.DPT'
filepathNormalizedN50 = '/Users/lukashein/Desktop/F-Praktikum-2/SolidStateOptics/RawData/StoN/StoN_res4_N50_normalized.DPT'
filepathNormalizedN75 = '/Users/lukashein/Desktop/F-Praktikum-2/SolidStateOptics/RawData/StoN/StoN_res4_N75_normalized.DPT'
filepathNormalizedN100 = '/Users/lukashein/Desktop/F-Praktikum-2/SolidStateOptics/RawData/StoN/StoN_res4_N100_normalized.DPT'

N10 = read_dpt_file(filepathNormalizedN10)
#N20 = read_dpt_file(filepathNormalizedN20)
#N50 = read_dpt_file(filepathNormalizedN50)
#N75 = read_dpt_file(filepathNormalizedN75)
N100 = read_dpt_file(filepathNormalizedN100)

# Plot für alle Datensätze
plot_data(N10)
#plot_data(N20)
#plot_data(N50)
#plot_data(N75)
plot_data(N100)

save_and_open()