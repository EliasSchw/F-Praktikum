from matplotlib import pyplot as plt
from DataPlotter import plot_data, save_and_open
from DataReader import read_dpt_file
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from Frauen import bügeln
from scipy.signal import find_peaks  # Added for peak detection

transmissionGaAsDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_doped_res03_N50_new_normalized.DPT')
#2000-9000
transmissionGaAsUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaAs_undoped_res03_N50_new_normalized.DPT')
#2000-10500
transmissionGaSbDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\GaSb_doped_res03_N50_new_normalized.DPT')
#GaSb nicht möglich wegen rauschen
transmissionSiUnDo = read_dpt_file(r'.\SolidStateOptics\RawData\Transmission_ex3\Si_undoped_res03_N50_normalized.DPT')
#SiUnDo 2000-8500

samplesOhneSiUn = [(transmissionSiUnDo, "Si Undoped", 530*10**-6, 0.0035, 2000, 8500),
        #(transmissionGaSbDo, "GaSb Doped", 500*10**-6, 69, 69, 69),#geht nicht
        (transmissionGaAsUnDo, "GaAs Undoped", 470*10**-6, 0.004, 2000, 10500),
        (transmissionGaAsDo, "GaAs Doped", 440*10**-6, 0.004, 2000, 9000)
    ]




def k_chopper(data, k_min, k_max):
    return [point for point in data if k_min <= point[0] <= k_max]

def plot_peaks(data, d, x_center, x_range, prominence=0.004, filename=""):
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]

    # Find peaks in the y-values
    peaks, _ = find_peaks(y_values, prominence=prominence, width=4)

    x_min = x_center - x_range / 2
    x_max = x_center + x_range / 2
    if x_center is not None and x_range is not None:
        peaks = [peak for peak in peaks if x_min <= x_values[peak] <= x_max]
    
    plt.plot([x for x in x_values if x_min <= x <= x_max], [y for x, y in zip(x_values, y_values) if x_min <= x <= x_max],
             '-', markersize=1)
    

    for i in range(len(peaks)):
       plt.plot(x_values[peaks[i]],y_values[peaks[i]], 'ro')
       #plt.plot(x_values[peaks[i]],0.55, 'ro')
       
       
    plt.xlabel(plt.gca().get_xlabel(), fontsize=15)
    plt.ylabel(plt.gca().get_xlabel(), fontsize=15)
    plt.tick_params(axis='both', labelsize=14)
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.tick_params(axis='both', length=6, width=1.2)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=15)
    plt.ylabel(r'Reflectivity R', fontsize=15)
    plt.xlim(left=x_min, right=x_max)
    plt.legend(fontsize=15) 
    save_and_open(filename=filename)

def calculate_n(data, d, x_center, windowSize, prominence=0.004):
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]

    peaks, _ = find_peaks(y_values, prominence=prominence, width=3)
    x_peaks = [x_values[peak] for peak in peaks]
    
    x_min = x_center - windowSize / 2
    x_max = x_center + windowSize / 2
    peaksInWindow = [peak for peak in peaks if x_min <= x_values[peak] <= x_max]
    
    periods = [- x_values[peaksInWindow[i + 1]] + x_values[peaksInWindow[i]] for i in range(len(peaksInWindow)-1)]
    period = sum(periods) / len(periods)
    print(f"nr. period {len(periods)} bei nu {x_center}")      
    return 1/(2*d*period*100)


def plot_the_ns(windowsize=200):
    colors = ['blue', 'green', 'red', 'purple', 'orange']  # Add a list of colors
    for idx, (data, label, d, prominence, nu_min, nu_max) in enumerate(samplesOhneSiUn):
        ns = []
        for nu in range(nu_min, nu_max, 500):
            ns.append(calculate_n(data, d, nu, windowsize, prominence=prominence))
            error = 0.02 * calculate_n(data, d, nu, windowsize, prominence=prominence)
            # ist die dicke der probe
            plt.errorbar(nu, ns[-1], yerr=error, fmt='o', color=colors[idx % len(colors)], capsize=5)  # Use color from the list
        plt.plot([nu for nu in range(nu_min, nu_max, 500)], ns, '.', label=label, markersize=10, color=colors[idx % len(colors)])  # Use color from the list

    plt.grid(True)
    plt.xlabel(plt.gca().get_xlabel(), fontsize=15)
    plt.ylabel(plt.gca().get_xlabel(), fontsize=15)
    plt.tick_params(axis='both', labelsize=14)
    plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    plt.tick_params(axis='both', length=6, width=1.2)
    plt.xlabel(r'wave number $\nu$ / ' + r'$cm^{-1}$', fontsize=15)
    plt.ylabel(r'refractive index n', fontsize=15)
    plt.legend(fontsize=15) 
    save_and_open("RefractiveIndicesDiscrete")


plot_peaks(transmissionSiUnDo, 500, 8000, 150, 0.004, filename="LewiPeak")
#plot_the_ns()


# Example usage
#average_period = calculate_n(transmissionGaAsDo, 500, x_center=2010, x_range=500)
#print(f"Average period {average_period}")

