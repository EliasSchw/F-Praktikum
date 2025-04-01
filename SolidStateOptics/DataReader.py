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