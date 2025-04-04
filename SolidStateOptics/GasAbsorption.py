from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
from scipy.signal import find_peaks
import numpy as np

# Dateipfade definieren
filepathNormalized = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_normalized.DPT'
filepathSample = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_sample.DPT'
filepathReference = './SolidStateOptics/RawData/gasAbsorption/gas_absorption_res03_N20_reference.DPT'

# Daten laden
dataSample = read_dpt_file(filepathSample)
dataReference = read_dpt_file(filepathReference)

#zum Überprüfen, ob die normalized Daten gleich meiner berechnnung sind
#divided_data=[]
#for s, r in zip(dataSample, dataReference):
#    divided_data.append([s[0], s[1] / r[1]])
#plot_data(divided_data)

# Daten laden und plotten
dataNormalized = read_dpt_file(filepathNormalized)
wave_numbers = [point[0] for point in dataNormalized]
transmissions = [point[1] for point in dataNormalized]
plot_data(dataNormalized)

# Manuelle Positionen der Dips und y-Werte für die Nummerierung
dip_positions = [2350, 3750, 5320, 7180]
y_positions = [0.64, 0.7, 0.965, 1]  # Angepasste y-Werte für die Nummerierung

for i, (dip, y_pos) in enumerate(zip(dip_positions, y_positions)):
    # Zahl direkt an den gegebenen Koordinaten platzieren
    plt.text(dip, y_pos, f'({i+1})', color='black', fontsize=10, ha='center')  # Nummerierung zentriert an den Koordinaten

plt.xlabel('wave number k / ' + r'$cm^{-1}$')
plt.ylabel('Transmission T')

save_and_open("GasAbsorption")