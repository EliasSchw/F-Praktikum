# c:\Users\schwa\OneDrive\EliasOneDrive\Uni\7. Semester\F-Praktikum\F-Praktikum\F-Praktikum\data_reader.py

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt


def plot_data(data):
    """
    Plots the data using matplotlib.

    Parameters:
        data (list): The data to plot, assumed to be a list of tuples or a 2D array.
    """
    # Extract x and y values from the data
    x_values = [point[0] for point in data]
    y_values = [point[1] for point in data]

    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, ',', markersize=100)  # Plot as single points
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plot of Data')
    plt.grid(True)
    plt.savefig('foo.png', dpi=400)
    from PIL import Image
    Image.open("foo.png").show()
    plt.clf()


def read_dpt_file(filepath):
    """
    Reads a .dpt file and returns its content as a list of tuples containing two float values.

    Args:
        filepath (str): Path to the .dpt file.

    Returns:
        list: A list of tuples, each containing two float values.
    """
    values = []
    try:
        with open(filepath, 'r') as file:
            for line in file:
                # Split the line into two values
                value1, value2 = map(float, line.split())
                values.append((value1, value2))
        return values
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# Example usage
filepath = r'C:/Users/schwa/OneDrive/EliasOneDrive/Uni/7. Semester/F-Praktikum/F-Praktikum/F-Praktikum/SolidStateOptics/RawData/ersterDatenCheck/StoN_res4_N10_normalized.DPT'
data = read_dpt_file(filepath)
plot_data(data)
#print(parsed_values)