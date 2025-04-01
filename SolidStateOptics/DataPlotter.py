# c:\Users\schwa\OneDrive\EliasOneDrive\Uni\7. Semester\F-Praktikum\F-Praktikum\F-Praktikum\data_reader.py

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from DataReader import read_dpt_file

plt.figure(figsize=(10, 6))

def plot_data(data):
    """
    Plots the data using matplotlib.

    Parameters:
        data (list): The data to plot, assumed to be a list of tuples or a 2D array.
    """
    # Extract x and y values from the data
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]
    

   
    plt.plot(x_values, y_values, '-', linewidth=0.3)  # Plot as single points
    plt.xlabel('wave number k / ' + r'$cm^{-1}$')
    plt.ylabel('Transmission T')
    plt.grid(True)
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.xlim(left=min(x_values), right=max(x_values))


def save_and_open(title="foo", filename="foo"):
    plt.title(title)
    plt.savefig(filename + '.png', dpi=600)
    from PIL import Image
    Image.open(filename + ".png").show()
    plt.clf()


filepathbeginning = r'C:/Users/schwa/OneDrive/EliasOneDrive/Uni/7. Semester/F-Praktikum/F-Praktikum/F-Praktikum/SolidStateOptics/RawData/'


file1 = r'C:\Users\schwa\OneDrive\EliasOneDrive\Uni\7. Semester\F-Praktikum\F-Praktikum\F-Praktikum\SolidStateOptics\RawData\Reflection_ex4\refl_GaAs_doped_res4_N50_normalized.DPT'
file3 = 'StoN_res4_N10_reference.DPT'

# Example usage
plot_data(read_dpt_file(file1))
#plot_data(read_dpt_file(filepathbeginning + file2))
#plot_data(read_dpt_file(filepathbeginning + file3))
save_and_open()