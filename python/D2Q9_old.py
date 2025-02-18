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
    """
    Class representation of the D2Q9 Lattice for the Lattice Boltzmann Method.
    """
    DIRECTIONS = 9  # Number of velocity directions

    # Velocity vectors for the D2Q9 lattice
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

    def __init__(self):
        """
        Initialize the D2Q9 Lattice object, including weights.
        """
        self.weights = self.compute_weights()
        self.velocities = self.compute_lattice_velocities()

    @staticmethod
    def compute_weights():
        """
        Compute the weights for the D2Q9 lattice.

        The weights are predefined constants that represent the proportion of particles 
        traveling in each direction within the D2Q9 model. These weights ensure 
        accurate simulation of fluid properties and preserve conservation laws.

        The directions and their corresponding weights are:
        - Stationary (center node): weight = 4/9
        - Horizontal and vertical directions: weight = 1/9 each
        - Diagonal directions: weight = 1/36 each

        Returns:
            numpy.ndarray: Array of weights for the D2Q9 lattice directions.
        """
        return np.array([
            4 / 9,  # Stationary
            *([1 / 9] * 4),  # Horizontal and vertical directions
            *([1 / 36] * 4)  # Diagonal directions
        ])

    def initialize_single_lattice_unit(self):
        """
        Initialize a single lattice unit for the D2Q9 model.

        This function sets up the initial values for a single unit in the lattice. The lattice 
        unit is represented as an array of 9 values (one for each velocity direction). 
        Each value corresponds to the particle distribution in that direction.

        Steps:
        1. Create an empty array of zeros for the 9 velocity directions.
        2. Assume an initial uniform density (rho = 1.0).
        3. Assign the equilibrium distribution for each direction based on the lattice weights.

        Returns:
            numpy.ndarray: A (1, 9) array representing the initial particle distribution 
            functions for a single lattice unit, rounded to three decimal places.
        """
        # Initialize an array for the distribution functions, one for each direction
        f = np.zeros((1, self.DIRECTIONS))

        # Set initial density to 1.0 (uniform distribution)
        rho = 1.0  # Initial density of the fluid at this lattice unit

        # Compute the initial distribution for each direction using the lattice weights
        for i in range(self.DIRECTIONS):
            f[0, i] = rho * self.weights[i]

        # Return the distribution functions rounded to 3 decimal places for clarity
        return np.round(f, 3)

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

        # Define velocities based on the equations
        velocities[0] = [0, 0]  # Stationary
        for i in range(1, 5):
            velocities[i] = [
                np.cos((i - 1) * np.pi / 2),
                np.sin((i - 1) * np.pi / 2)
            ]
        for i in range(5, 9):
            velocities[i] = np.sqrt(2) * np.array([
                np.cos((i - 5) * np.pi / 2 + np.pi / 4),
                np.sin((i - 5) * np.pi / 2 + np.pi / 4)
            ])

        return np.round(velocities, 3)


    def equilibrium_distribution_single_unit(self, rho, ux, uy):
        """
        Compute the equilibrium distribution for a single lattice unit in the D2Q9 model.

        This function calculates the equilibrium distribution, which represents the idealized
        particle distribution in each direction based on the local fluid density and velocity.
        The equilibrium distribution ensures consistency with the macroscopic fluid equations.

        Steps:
        1. Initialize an array of zeros for the distribution values in 9 directions.
        2. For each direction, compute the dot product of the velocity vector with the
           macroscopic fluid velocity (ux, uy).
        3. Use the LBM equilibrium formula to calculate the distribution for each direction:
           f_eq[i] = w[i] * rho * (1 + (xi[i] . u) / c_s^2 + (xi[i] . u)^2 / (2*c_s^4) - u.u / (2*c_s^2))
        4. Return the computed distribution rounded to 3 decimal places.

        Args:
            rho (float): Local fluid density.
            ux (float): Velocity in the x-direction.
            uy (float): Velocity in the y-direction.

        Returns:
            numpy.ndarray: Equilibrium distribution function for a single lattice unit.
        """
        # Initialize an array to store the equilibrium distribution for 9 directions
        feq = np.zeros(self.DIRECTIONS)

        # Speed of sound squared
        c_s2 = 1 / 3

        # Loop through each velocity direction
        for i, xi in enumerate(self.velocities):
            # Compute the dot product of the velocity vector and the macroscopic velocity
            cu = (xi[0] * ux + xi[1] * uy)

            # Normalize the dot product by the speed of sound squared
            cu = cu / c_s2

            # Calculate the equilibrium distribution for this direction using the LBM formula
            feq[i] = self.weights[i] * rho * (
                1 + cu + 0.5 * cu**2 - 0.5 * (ux**2 + uy**2) / c_s2
            )

        # Return the computed equilibrium distribution rounded for readability
        return np.round(feq, 3)

    def compute_density(self, f):
        """
        Compute the macroscopic density from the particle distribution functions.

        This function calculates the macroscopic density at a lattice unit based on the
        particle distribution functions. The density is the sum of the distribution values
        in all directions.

        Args:
            f (numpy.ndarray): Array of particle distribution functions for a lattice unit.

        Returns:
            float: Macroscopic density computed from the distribution functions.
        """
        return np.sum(f)
    
    def compute_momentum(self, xi, f,):
        """
        Compute the momentum of a single lattice unit.

        Parameters:
        xi (numpy.ndarray): The velocity vectors.
        f (numpy.ndarray): The distribution functions.

        Returns:
        numpy.ndarray: The computed momentum as a dot product of xi and f.
        """
        return np.dot(xi, f)

