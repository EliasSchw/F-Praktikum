from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file


file1 = r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_doped_res03_N50_new_normalized.DPT'

ausschnitt = []
for s in read_dpt_file(file1):
    if s[0] > 2000 and s[0] < 2200:
        ausschnitt.append([s[0], s[1]])

plot_data(ausschnitt)
save_and_open()