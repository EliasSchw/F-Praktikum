from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
from Glätteisen import bügeln


reflection = read_dpt_file(r'.\SolidStateOptics\RawData\Reflection_ex4\refl_GaAs_doped_res4_N50_normalized.DPT')
transmission = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex4\GaAs_doped_res4_N50_normalized.DPT')



# ausschnitt = []
# for s in read_dpt_file(file1):
#     if s[0] > 2000 and s[0] < 2200:
#         ausschnitt.append([s[0], s[1]])


plot_data(bügeln(reflection,1))
plot_data(bügeln(reflection,0.025))
plt.xlabel('wave number k / ' + r'$cm^{-1}$')
plt.ylabel('Reflection R')
save_and_open()

plot_data(transmission)
plt.xlabel('wave number k / ' + r'$cm^{-1}$')
plt.ylabel('Transmission T')
save_and_open()