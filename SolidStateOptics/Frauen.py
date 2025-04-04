def bügeln(data, alpha = 0.5):
    """
    Applies exponential smoothing to the given data. 
    Higher alpha values give more weight to recent observations.
    Lower alpha values give more weight to older observations.

    Args:
        data (list): A list of tuples, where each tuple contains two float values (x, y).
        alpha (float): Smoothing factor, where 0 < alpha <= 1.

    Returns:
        list: A list of tuples with smoothed y-values.
    """
    if not 0 < alpha <= 1:
        raise ValueError("Alpha must be between 0 and 1.")

    smoothed_data = []
    if data:
        # Initialize the first smoothed value with the average value of the first 10 y-values
        initial_values = [y for _, y in data[:100]]
        smoothed_y = sum(initial_values) / len(initial_values)
        smoothed_data.append((data[0][0], smoothed_y))

        # Apply exponential smoothing to the rest of the data
        for x, y in data[1:]:  #alle bis auf den ersten
            smoothed_y = alpha * y + (1 - alpha) * smoothed_y
            smoothed_data.append((x, smoothed_y))

    return smoothed_data

def average_filter(data, window_size=3):
    """
    Applies a moving average filter to smooth noise in the given data, considering edge cases.

    Args:
        data (list): A list of tuples, where each tuple contains two float values (x, y).
        window_size (int): The size of the moving window for averaging.

    Returns:
        list: A list of tuples with smoothed y-values.
    """
    if window_size < 1:
        raise ValueError("Window size must be at least 1.")

    smoothed_data = []
    half_window = window_size // 2

    for i in range(len(data)):
        # Determine the start and end indices of the window
        start_idx = max(0, i - half_window)
        end_idx = min(len(data), i + half_window + 1)

        # Handle edges by extending the window with the nearest values
        window = [data[j][1] if 0 <= j < len(data) else data[i][1] for j in range(start_idx, end_idx)]
        avg_y = sum(window) / len(window)

        # Append the smoothed value
        smoothed_data.append((data[i][0], avg_y))

    return smoothed_data