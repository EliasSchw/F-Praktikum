from SolidStateOptics.DataPlotter import plot_data, save_and_open
from SolidStateOptics.DataReader import read_dpt_file

filepathbeginning = './SolidStateOptics/RawData/'
file1 = 'Reflection_ex4/refl_GaAs_doped_res4_N50_normalized.DPT'

plot_data(read_dpt_file(file1))
save_and_open()