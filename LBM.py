import numpy as np

def compute_equilibrium(rho, u):
    """
    Compute the equilibrium distribution function (f_eq).

    Args:
        rho (float): Macroscopic density.
        u (numpy.ndarray): Macroscopic velocity vector.

    Returns:
        numpy.ndarray: Equilibrium distribution functions for all directions.
    """

    ksi = np.array([
        [0, 0],  # Stationary
        [1, 0],  # Right
        [0, 1],  # Up
        [-1, 0], # Left
        [0, -1], # Down
        [1, 1],  # Top-right
        [-1, 1], # Top-left
        [-1, -1],# Bottom-left
        [1, -1]  # Bottom-right
    ])

    weights = np.array([
            4 / 9,  # Stationary
            *([1 / 9] * 4),  # Horizontal and vertical directions
            *([1 / 36] * 4) ])

    cs2 = 1 / 3  # Square of speed of sound

    feq = np.zeros(9)
    for i, xi in enumerate(ksi):
        cu = np.dot(xi, u)  # Dot product of velocity vector and u
        feq[i] = weights[i] * rho * (
            1 + cu / cs2 + 0.5 * (cu**2) / cs2**2 - 0.5 * (np.dot(u, u)) / cs2
        )
    return np.round(feq, 3)