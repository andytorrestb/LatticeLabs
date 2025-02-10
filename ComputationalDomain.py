import numpy as np

class ComputationalDomain:
    
    def __init__(self):
        """
        Initialize the ComputationalDomain object.
        """
        pass

class PoiseuilleFlow(ComputationalDomain):

    def __init__(self, nx, ny, dx, dy):
        """
        Initialize the PoiseuilleFlow object.

        Args:
            nx (int): Number of lattice nodes in the x-direction.
            ny (int): Number of lattice nodes in the y-direction.
            dx (float): Spacing between lattice nodes in the x-direction.
            dy (float): Spacing between lattice nodes in the y-direction.
        """
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        # self.initialize()

    def set_lattice(self, lattice):
        """
        Set the lattice object for the computational domain.

        Args:
            lattice (Lattice): Lattice object (e.g., D2Q9).
        """
        self.lattice = lattice

    def set_inlet_density(self, rho_in):
        """
        Set the inlet density for the computational domain.

        Args:
            rho_in (float): Inlet density value.
        """
        self.rho_in = rho_in 

    def set_outlet_zerogradient(self):
        """
        Set the boundary conditions for the computational domain.

        Args:
            boundary_conditions (BoundaryConditions): BoundaryConditions object.
        """
        self.rho_out = 0.0

    def initialize(self):
        """
        Initialize the Poiseuille flow computational domain.
        """

        assert(hasattr(self, 'lattice')), "Lattice object not set."


        directions = self.lattice.DIRECTIONS

        # Initialize the f and f_eq arrays
        self.f = np.zeros((self.ny, self.nx, directions))
        self.f_eq = np.zeros((self.ny, self.nx, directions))

        # Initialize the density and velocity fields
        self.rho = np.ones((self.ny, self.nx)) * self.rho_in
        self.u = np.zeros((self.ny, self.nx, 2))


        # Initialize the equilibrium distribution
        for j in range(self.ny):
            for i in range(self.nx):
                self.f_eq[j, i] = self.lattice.compute_equilibrium(self.rho[j, i], self.u[j, i])
                # print(self.f_eq[j, i])
