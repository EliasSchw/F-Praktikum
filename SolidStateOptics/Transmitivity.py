from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from Frauen import bügeln


transmission = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_undoped_res03_N50_new_normalized.DPT')

x1_inset = 6000
x2_inset = 6100



def k_chopper(data, k_min, k_max):
    return [point for point in data if k_min <= point[0] <= k_max]


plot_data(transmission)

plt.axvline(x=x1_inset, color='red', linestyle='-', label='k=6000', linewidth=0.5)
plt.axvline(x=x2_inset, color='red', linestyle='-', label='k=6100', linewidth=0.5)

#plt.legend()

# Create an inset plot

# Chop the data for the inset
chopped_data = k_chopper(transmission, x1_inset, x2_inset)

# Create the inset axes
ax = plt.gca()
ax_inset = ax.inset_axes([0.15, 0.08, 0.50, 0.50])
ax_inset.plot([point[0] for point in chopped_data], [point[1] for point in chopped_data], color='blue')
#ax_inset.set_title("Inset")
#ax_inset.set_xlabel(r'$\nu$')
#ax_inset.set_ylabel('T')
ax_inset.tick_params(axis='both', direction='in')
ax_inset.tick_params(axis='both', which='both', direction='in', top=True, right=True)



plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=18)
plt.ylabel(r'transmission T', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
save_and_open(filename="foo")#Transmission_High_Res_GaAs_Doped_N50")