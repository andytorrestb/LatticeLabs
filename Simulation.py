import numpy as np

class Simulation():
    # The Simulation class is responsible for running the simulation.
    # It contains the main loop of the simulation.
    # The Simulation class is initialized with a D2Q9 object and a ComputationalDomain object.
    # The Simulation class has a method run that runs the simulation for a given number of time steps.

    def __init__(self, lattice, computational_domain):
        """
        Initialize the Simulation object.

        Args:
            d2q9 (D2Q9): D2Q9 lattice object.
            computational_domain (ComputationalDomain): Computational domain object.
        """
        self.lattice = lattice
        self.computational_domain = computational_domain
        return

    def update(self):
        """
        Update the simulation for one time step.
        """

        f_new = np.zeros((self.computational_domain.ny, self.computational_domain.nx, self.lattice.DIRECTIONS))
        f = self.computational_domain.f

        n_x = self.computational_domain.nx
        n_y = self.computational_domain.ny  

        x_max = n_x - 1
        y_min = n_y - 1

        for j in range(self.computational_domain.ny):
            for i in range(self.computational_domain.nx):

                ## Streaming
                match [j, i]:
                    case [0, 0]: # Top left corner
                        pass
                    case [0, x_max]: # Top right corner
                        pass
                    case [0, _]: # Interior top surface
                        pass
                    case [y_min, 0]: # Bottom left corner
                        pass
                    case [y_min, x_max]: # Bottom right corner
                        pass
                    case [y_min, _]: # Interior bottom surface
                        pass
                    case [_, 0]: # Left surface
                        pass
                    case [_, x_max]: # Right surface
                        pass
                    case _: # Interior nodes
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 2] = f[j+1, i, 2]
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 4] = f[j-1, i, 4]
                        f_new[j, i, 5] = f[j+1, i-1, 5]
                        f_new[j, i, 6] = f[j+1, i+1, 6]
                        f_new[j, i, 7] = f[j-1, i+1, 7]
                        f_new[j, i, 8] = f[j-1, i-1, 8]

                        # Collision
                        rho = self.computational_domain.rho[j, i]
                        u = self.computational_domain.u[j, i]
                        self.computational_domain.f_eq[j, i] = self.lattice.compute_equilibrium(rho, u)
                        self.computational_domain.f[j, i] = self.computational_domain.rho_in * self.lattice.weights
                        self.computational_domain.rho[j, i] = np.sum(self.computational_domain.f[j, i]) + 1

        # print("Updating simulation...")
        

    def run(self, time_steps):
        """
        Run the simulation for a given number of time steps.

        Args:
            time_steps (int): Number of time steps to run the simulation.
        """

        results = []

        u = self.computational_domain.u
        rho = self.computational_domain.rho
        for t in range(time_steps):
            # self.update()
            rho = rho + 1
            result = {"time": t, "density": rho, "velocity": u}
            results.append(result)
        return results
