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
        # Initialize the first smoothed value with the first y-value
        smoothed_y = data[0][1]
        smoothed_data.append((data[0][0], smoothed_y))

        # Apply exponential smoothing to the rest of the data
        for x, y in data[1:]:  #alle bis auf den ersten
            smoothed_y = alpha * y + (1 - alpha) * smoothed_y
            smoothed_data.append((x, smoothed_y))

    return smoothed_data
