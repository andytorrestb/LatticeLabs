import numpy as np
from LBM import compute_equilibrium

# Define computational domain parameters
n_x = 100 
n_y = int(n_x / 2)
dx = 1.0
dy = 1.0

rho_in = 5.0

# Define LBM parameters

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

tau = 1.0

# Initialize computational domain
f_eq = np.zeros((n_y, n_x, 9))
f = np.zeros((n_y, n_x, 9))
rho = np.ones((n_y, n_x)) * rho_in
u = np.ones((n_y, n_x, 2)) * rho_in
u_x = u[:, :, 0]
u_y = u[:, :, 1]

for j in range(n_y):
    for i in range(n_x):
        f_eq[j, i] = compute_equilibrium(rho[j, i], u[j, i])
        for k in range(9):
            f[j, i, k] = rho_in * weights[k]

# print(f_eq)
# print(f)
# Run the simulation
f_new = np.zeros((n_y, n_x, 9))
x_max = n_x - 1
y_min = n_y - 1

T=2000

for t in range(T):
    ## Streaming
    for j in range(n_y):
        for i in range(n_x):

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
        ## Collision
        f = f_new - (f_new - f_eq) / tau

        # Update macroscopic variables
        rho = np.sum(f, axis=2)
        print("ksi x", ksi[:, 0])
        print("ksi y", ksi[:, 1])


        u_x = np.dot(f, ksi[:, 0]) / rho
        u_y = np.dot(f, ksi[:, 1]) / rho
        print("f", f)
        print("u_x", u_x)
        print("u_y", u_y)
# print(u_y)

# Visualization
import matplotlib.pyplot as plt

# Visualization
X, Y = np.meshgrid(np.arange(n_x), np.arange(n_y))
plt.quiver(X, Y, u_x, u_y)
plt.title('Velocity field')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('velocity_plot.png')

