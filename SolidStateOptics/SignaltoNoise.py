from matplotlib import pyplot as plt
from DataReader import read_dpt_file
import numpy as np
from scipy.stats import linregress
import os, sys
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-2]))
from macroswriter import writeLatexMacro

filepathNormalizedN10 = './SolidStateOptics/RawData/StoN/StoN_res4_N10_normalized.DPT'
filepathNormalizedN20 = './SolidStateOptics/RawData/StoN/StoN_res4_N20_normalized.DPT'
filepathNormalizedN50 = './SolidStateOptics/RawData/StoN/StoN_res4_N50_normalized.DPT'
filepathNormalizedN75 = './SolidStateOptics/RawData/StoN/StoN_res4_N75_normalized.DPT'
filepathNormalizedN100 = './SolidStateOptics/RawData/StoN/StoN_res4_N100_normalized.DPT'
#plotrange
yMin = 0.98
yMax = 1.02
yMinInset = 0.998
yMaxInset = 1.002
plt.figure(figsize=(10, 6))
def plot_data(data, label=''):
    """
    Plots the data using matplotlib.

    Parameters:
        data (list): The data to plot, assumed to be a list of tuples or a 2D array.
    """
    # Extract x and y values from the data
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]
    
    plt.plot(x_values, y_values, '-', linewidth=0.9, label=label)  # Plot as single points
    plt.grid(True)
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.xlim(left=min(x_values), right=max(x_values))
    plt.ylim(yMin, yMax)

def save_and_open(filename="SignalToNoise", title=""):
    plt.legend(fontsize=12)  # Legende vergrößert
    plt.xlabel(r'wave number $cm^{-1}$')  # x-Achse im LaTeX-Stil benannt
    plt.ylabel(r'transmission $T$')  # y-Achse im LaTeX-Stil benannt

    # Inset-Plot hinzufügen
    ax = plt.gca()
    inset_ax = ax.inset_axes([0.15, 0.05, 0.35, 0.35])  # Position und Größe des Insets [x, y, Breite, Höhe]

    # Daten für das Inset plotten
    for line in ax.get_lines():
        x_data = line.get_xdata()
        y_data = line.get_ydata()

        # Bereich für das Inset auswählen
        mask = (x_data >= 6000) & (x_data <= 7000)
        x_inset = x_data[mask]
        y_inset = y_data[mask]

        # Mittelwert der beiden Signale berechnen und verschieben
        y_mean = np.mean(y_inset)
        y_shifted = y_inset - y_mean + 1.0  # Verschieben, sodass der Mittelwert bei 1 liegt

        inset_ax.plot(x_inset, y_shifted, line.get_linestyle(), linewidth=0.9, label=line.get_label())

    # Rechte y-Achse für das Inset hinzufügen
    inset_ax_right = inset_ax.twinx()
    inset_ax_right.set_ylim(inset_ax.get_ylim())  # Gleiche Skalierung wie die linke Achse
    inset_ax_right.tick_params(axis='y', direction='in', which='both', right=True, labelright=True, labelleft=False)

    # Einheitliche Achsenbeschriftung und Skalierung
    inset_ax.set_xlim(6000, 7000)
    inset_ax.set_ylim(0.998, 1.002)  # Einheitliche y-Achsenbegrenzung
    inset_ax.set_ylabel(r'transmission $T$', fontsize=10)  # Beschriftung für die linke y-Achse
    inset_ax_right.set_ylabel(r'transmission $T$', fontsize=10)  # Gleiche Beschriftung für die rechte y-Achse

    plt.savefig('.\\Paper\\Images\\'+filename + '.png', dpi=600)
    from PIL import Image
    Image.open(".\\Paper\\Images\\"+filename + ".png").show()
    plt.clf()

N10 = read_dpt_file(filepathNormalizedN10) 
N20 = read_dpt_file(filepathNormalizedN20)
N50 = read_dpt_file(filepathNormalizedN50)
N75 = read_dpt_file(filepathNormalizedN75)
N100 = read_dpt_file(filepathNormalizedN100)


datasets = [(N10, 'N=10'), (N20, 'N=20'), (N50, 'N=50'), (N75, 'N=75'), (N100, 'N=100')]
x_min = 6000
x_max = 7000

def calculateSNRWithError(datasets, x_min, x_max):
    """
    Berechnet die Signal-to-Noise-Ratio (SNR) und die zugehörigen Fehler für eine Liste von Datensätzen.

    Parameters:
        datasets (list): Liste von Tupeln, die die Daten und Labels enthalten.
        x_min (float): Untere Grenze des x-Bereichs.
        x_max (float): Obere Grenze des x-Bereichs.

    Returns:
        tuple: Liste von SNR-Werten, zugehörigen Fehlern und Labels.
    """
    snr_results = []
    snr_errors = []
    snr_labels = []
    for data, label in datasets:
        x_values = [row[0] for row in data]  # Extrahiere die x-Werte
        y_values = [row[1] for row in data]  # Extrahiere die y-Werte
        
        # Werte im Bereich filtern
        filtered_y = [y for x, y in zip(x_values, y_values) if x_min < x < x_max]
        
        # Mittelwert und Standardabweichung berechnen
        mean_value = np.mean(filtered_y)
        std_dev = np.std(filtered_y)
        
        # SNR berechnen
        snr = mean_value / std_dev if std_dev != 0 else float('inf')
        snr_results.append(snr)
        snr_labels.append(label)
        
        # Gaußsche Fehlerfortpflanzung für den Fehler des SNR
        error = snr * np.sqrt((std_dev / mean_value)**2 + (1 / len(filtered_y)))
        snr_errors.append(error)
    
    return snr_results, snr_errors, snr_labels

def plotSNR(datasets, x_min, x_max):
    """
    Plottet die Signal-to-Noise-Ratio (SNR) mit Fehlerbalken und führt eine lineare Regression im Log-Log-Raum durch.

    Parameters:
        datasets (list): Liste von Tupeln, die die Daten und Labels enthalten.
        x_min (float): Untere Grenze des x-Bereichs.
        x_max (float): Obere Grenze des x-Bereichs.
    """
    snr_results, snr_errors, snr_labels = calculateSNRWithError(datasets, x_min, x_max)
    
    # Konvertiere Labels in numerische Werte (z. B. 10, 20, 50, 75, 100)
    x_numeric = np.array([int(label.split('=')[1]) for _, label in datasets])
    y_numeric = np.array(snr_results)
    y_errors = np.array(snr_errors)
    
    # Log-Log-Transformation
    log_x = np.log(x_numeric)
    log_y = np.log(y_numeric)
    
    # Lineare Regression im Log-Log-Raum mit linregress
    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)
    regression_line = slope * log_x + intercept

    # Plot
    plt.figure(figsize=(8, 5))
    plt.errorbar(x_numeric, y_numeric, yerr=y_errors, fmt='o', color='blue', capsize=2.5, label='SNR Data Points')
    plt.plot(x_numeric, np.exp(regression_line), color='red', linestyle='--', label=f'linear regression')
    plt.xscale('log')  # Logarithmische Skalierung der x-Achse
    plt.yscale('log')  # Logarithmische Skalierung der y-Achse
    plt.xlabel('Number of Scans (N)')
    plt.ylabel('Signal-to-Noise Ratio (SNR)')
    #plt.title('Log-Log Plot of Signal-to-Noise Ratio with Errors')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)  # Ticks an allen Seiten nach innen
    plt.legend()

    # Exponenten und Fehler ausgeben und in LaTeX-Makro schreiben
    print(f"Exponent (x) from Log-Log Fit: {slope}")
    print(f"Standard Error of Exponent: {std_err}")
    writeLatexMacro('stn', slope, std_err)
    plt.savefig('.\\Paper\\Images\\'+'STNLog' + '.png', dpi=600)
    from PIL import Image
    Image.open(".\\Paper\\Images\\"+'STNLog' + ".png").show()
    plt.clf()


# Plot für alle Datensätze
plot_data(N10)
#plot_data(N20)
#plot_data(N50)
#plot_data(N75)
plot_data(N100)
save_and_open()


plotSNR(datasets, x_min, x_max)
#save_and_open()
