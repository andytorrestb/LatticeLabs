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
        rho_in = self.computational_domain.rho_in

        x_max = n_x - 1
        y_min = n_y - 1

        for j in range(self.computational_domain.ny):
            for i in range(self.computational_domain.nx):

                ## Streaming
                match [j, i]:
                    case [0, 0]: # Top left corner

                        # Known PDF values.
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 2] = f[j+1, i, 2]
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 6] = f[j+1, i+1, 6]

                        # Hard coded solutions for Unknown PDF values.
                        f_new[j, i, 1] = f[j, i, 3]
                        f_new[j, i, 4] = f[j, i, 2]
                        f_new[j, i, 5] = 0.5*(rho_in - f[j, i, 0]) - f[j, i, 2] - f[j, i, 3] - f[j, i, 6]
                        f_new[j, i, 7] = f[j, i, 5]
                        f_new[j, i, 8] = f[j, i, 6]

                    case [0, x_max]: # Top right corner
                         # Apply zero gradient boundary condition.
                        f_new[j, i, 0] = f[j, i-1, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 2] = f[j, i-1, 2]
                        f_new[j, i, 3] = f[j, i-1, 3]
                        f_new[j, i, 4] = f[j, i-1, 4]
                        f_new[j, i, 5] = f[j, i-1, 5]
                        f_new[j, i, 6] = f[j, i-1, 6]
                        f_new[j, i, 7] = f[j, i-1, 7]
                        f_new[j, i, 8] = f[j, i-1, 8]
                    case [0, _]: # Interior top surface
                        # Known PDF values.
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 2] = f[j+1, i, 2]
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 5] = f[j+1, i-1, 5]
                        f_new[j, i, 6] = f[j+1, i+1, 6]

                        # Hard coded solutions for Unknown PDF values.
                        f_new[j, i, 4] = f[j, i, 2]
                        f_new[j, i, 7] = 0.5*(f[j, i, 1] - f[j, i, 3]) + f[j, i, 5]
                        f_new[j, i, 8] = 0.5*(f[j, i, 3] - f[j, i, 1]) + f[j, i, 6]
                    case [y_min, 0]: # Bottom left corner
                        # Known PDF values.
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 4] = f[j-1, i, 4]
                        f_new[j, i, 7] = f[j-1, i+1, 7]

                        # Hard coded solutions for Unknown PDF values.
                        f_new[j, i, 1] = f[j, i, 3]
                        f_new[j, i, 2] = f[j, i, 4]
                        f_new[j, i, 5] = f[j, i, 7]
                        f_new[j, i, 6] = 0.5*(rho_in - f[j, i, 0]) - f[j, i, 3] - f[j, i, 4] - f[j, i, 7]
                        f_new[j, i, 8] = f[j, i, 6]
                        
                    case [y_min, x_max]: # Bottom right corner
                         # Apply zero gradient boundary condition.
                        f_new[j, i, 0] = f[j, i-1, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 2] = f[j, i-1, 2]
                        f_new[j, i, 3] = f[j, i-1, 3]
                        f_new[j, i, 4] = f[j, i-1, 4]
                        f_new[j, i, 5] = f[j, i-1, 5]
                        f_new[j, i, 6] = f[j, i-1, 6]
                        f_new[j, i, 7] = f[j, i-1, 7]
                        f_new[j, i, 8] = f[j, i-1, 8]
                        pass
                    case [y_min, _]: # Interior bottom surface
                        # Simply stream data from neighboring nodes.
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 4] = f[j-1, i, 4]
                        f_new[j, i, 7] = f[j-1, i+1, 7]
                        f_new[j, i, 8] = f[j-1, i-1, 8]

                        # Hard coded solutions for Unknown PDF values.
                        f_new[j, i, 2] = f[j, i, 4]
                        f_new[j, i, 5] = 0.5*(f[j, i, 3] - f[j, i, 1]) + f[j, i, 7]
                        f_new[j, i, 6] = 0.5*(f[j, i, 1] - f[j, i, 3]) + f[j, i, 8]
                    case [_, 0]: # Left surface
                        # Known PDF values.
                        f_new[j, i, 0] = f[j, i, 0] 
                        f_new[j, i, 2] = f[j+1, i, 2]
                        f_new[j, i, 3] = f[j, i+1, 3]
                        f_new[j, i, 4] = f[j-1, i, 4]
                        f_new[j, i, 6] = f[j+1, i+1, 6]
                        f_new[j, i, 7] = f[j-1, i+1, 7]

                        # Hard coded solutions for Unknown PDF values.
                        U_in = (f[j, i, 0] + f[j, i, 2] + f[j, i, 4] + 2*(f[j, i, 3] + f[j, i, 6] + f[j, i, 7])) / rho_in
                        f_new[j, i, 1] = f[j, i, 3] + (2/3) * rho_in * U_in
                        f_new[j, i, 5] = f[j, i, 7] + 0.5 * (f[j, i, 4] - f[j, i, 2]) + (1/6) * rho_in * U_in
                        f_new[j, i, 8] = f[j, i, 6] - 0.5 * (f[j, i, 4] - f[j, i, 2]) + (1/6) * rho_in * U_in
                    case [_, x_max]: # Right surface
                         # Apply zero gradient boundary condition.
                        f_new[j, i, 0] = f[j, i-1, 0] 
                        f_new[j, i, 1] = f[j, i-1, 1]
                        f_new[j, i, 2] = f[j, i-1, 2]
                        f_new[j, i, 3] = f[j, i-1, 3]
                        f_new[j, i, 4] = f[j, i-1, 4]
                        f_new[j, i, 5] = f[j, i-1, 5]
                        f_new[j, i, 6] = f[j, i-1, 6]
                        f_new[j, i, 7] = f[j, i-1, 7]
                        f_new[j, i, 8] = f[j, i-1, 8]
                    case _: # Interior nodes
                        # Simply stream data from neighboring nodes.
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
        self.computational_domain.f = f_new - (f_new - self.computational_domain.f_eq) / self.lattice.tau
    
        # Update the macroscopic density and velocity fields post-collision.

        ksi = self.lattice.LATTICE_VELOCITES
        for j in range(n_y):
            for i in range(n_x):
        
                self.computational_domain.rho[j, i] = np.sum(self.computational_domain.f[j, i])
                self.computational_domain.u_x[j, i] = np.sum(f[j, i]*ksi[:, 0])/self.computational_domain.rho[j, i]
                self.computational_domain.u_y[j, i] = np.sum(f[j, i]*ksi[:, 1])/self.computational_domain.rho[j, i]
                self.computational_domain.u[j, i] = np.sqrt(self.computational_domain.u_x[j, i]**2 + self.computational_domain.u_y[j, i]**2)

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
            self.update()
            rho = self.computational_domain.rho
            u = self.computational_domain.u
            result = {"time": t, "density": rho, "velocity": u}
            results.append(result)
        return results
