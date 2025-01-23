import numpy as np

class Lattice:
    """
    Base class for lattice models in the Lattice Boltzmann Method.
    """

    def __init__(self):
        """
        Initialize the Lattice object.
        """
        pass

class D2Q9(Lattice):

    # D2Q9 Lattice Parameters

    DIRECTIONS = 9  # Number of velocity directions

    # Directionality Velocity vectors for the D2Q9 lattice
    VELOCITIES = np.array([
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

    # Weights for the D2Q9 lattice
    WEIGHTS = np.array([
        4/9,  # Stationary
        1/9, 1/9, 1/9, 1/9,  # Cardinal directions
        1/36, 1/36, 1/36, 1/36  # Diagonal directions
    ])

    CS = 1 / np.sqrt(3)  # Speed of sound in lattice units

    def __init__(self):
        """
        Initialize the BGKCollision class with lattice parameters.

        Args:
            velocities (numpy.ndarray): Discrete velocity vectors (e.g., D2Q9).
            weights (numpy.ndarray): Weight for each lattice direction.
            cs (float): Speed of sound in lattice units.
        """
        self.velocities = self.VELOCITIES
        self.weights = self.WEIGHTS
        self.cs = self.CS
        self.cs2 = self.cs**2  # Square of speed of sound
        self.compute_lattice_velocities()

    def compute_density(self, f):
        """
        Compute the macroscopic density (rho).

        Args:
            f (numpy.ndarray): Distribution functions for a node.

        Returns:
            float: Macroscopic density.
        """
        return np.sum(f)

    def compute_lattice_velocities(self):
        """
        Compute the lattice velocities based on the D2Q9 model.

        The velocities are calculated according to the equations in the LBM framework for D2Q9:
        - The stationary velocity is [0, 0]
        - For directions 1 to 4, velocities are [cos((i-1) * pi / 2), sin((i-1) * pi / 2)]
        - For directions 5 to 8, velocities are sqrt(2) * [cos((i-5) * pi / 2 + pi / 4), sin((i-5) * pi / 2 + pi / 4)]

        Returns:
            numpy.ndarray: Array of lattice velocities for the D2Q9 model.
        """
        velocities = np.zeros((self.DIRECTIONS, 2))

        velocities[0] = [0, 0]  # Stationary
    
        for i in range(1, 5): # Cardinal directions
            velocities[i] = [
                np.cos((i - 1) * np.pi / 2),
                np.sin((i - 1) * np.pi / 2)
            ]

        for i in range(5, 9): # Diagonal directions
            velocities[i] = np.sqrt(2) * np.array([
                np.cos((i - 5) * np.pi / 2 + np.pi / 4),
                np.sin((i - 5) * np.pi / 2 + np.pi / 4)
            ])

        self.LATTICE_VELOCITES = np.round(velocities, 3)
        return np.round(velocities, 3)

    def compute_velocity(self, f, rho):
        """
        Compute the macroscopic velocity (u).

        Args:
            f (numpy.ndarray): Distribution functions for a node.
            rho (float): Macroscopic density.

        Returns:
            numpy.ndarray: Macroscopic velocity vector (u).
        """
        return np.dot(f, self.LATTICE_VELOCITES) / rho

    def compute_equilibrium(self, rho, u):
        """
        Compute the equilibrium distribution function (f_eq).

        Args:
            rho (float): Macroscopic density.
            u (numpy.ndarray): Macroscopic velocity vector.

        Returns:
            numpy.ndarray: Equilibrium distribution functions for all directions.
        """
        feq = np.zeros(len(self.velocities))
        for i, xi in enumerate(self.velocities):
            cu = np.dot(xi, u)  # Dot product of velocity vector and u
            feq[i] = self.weights[i] * rho * (
                1 + cu / self.cs2 + 0.5 * (cu**2) / self.cs2**2 - 0.5 * (np.dot(u, u)) / self.cs2
            )
        return np.round(feq, 3)

def main():

    # Table of particle distribution functions (PDF) for each node in the lattice
    nodes = {  
        "A": np.array([1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16]),
        "B": np.array([1.67, 0.42, 0.42, 0.42, 0.42, 0.1, 0.11, 0.1, 0.11]),
        "C": np.array([1.66, 0.5, 0.42, 0.35, 0.42, 0.12, 0.09, 0.08, 0.13])
    }

    # Initialize the D2Q9 lattice
    d2q9 = D2Q9()


    # Compute and print results for each node
    for node, f in nodes.items():
        rho = d2q9.compute_density(f)
        u = d2q9.compute_velocity(f, rho)
        feq = d2q9.compute_equilibrium(rho, u)

        print(f"Node {node}:")
        print(f"  Density (rho): {rho:.3f}")
        print(f"  Velocity (u): {u}")
        print(f"  Equilibrium Distribution (f_eq): {feq}")
        print()

if __name__ == "__main__":
    main()
