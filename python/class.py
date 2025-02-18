import numpy as np
import D2Q9
# Define parameters

# Computational Domain related
n_x = 100
n_y = int(n_x / 2)

dx = 1.0
dy = 1.0



# LBM related
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
lattice = D2Q9.D2Q9(tau)

# Initialize computational domain
f_eq = np.zeros((n_y, n_x, 9))
f = np.zeros((n_y, n_x, 9))

for j in range(n_y):
    for i in range(n_x):
        f_eq[j, i] = lattice.compute_equilibrium(rho[j, i], u_x[j, i], u_y[j, i])
        for k in range(9):

            f[k] = rho_in*weights[k]

f_new = np.zeros((n_y, n_x, 9))

rho_in = 5.0
rho = np.ones((n_y, n_x))*rho_in
u_x = np.zeros((n_y, n_x))
u_x = np.zeros((n_y, n_x))

# Solving 
for j in range(n_y):
    for i in range(n_x):

## Streaming
        if j == 0: # Top surface

            if i == 0: # Top left corner
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 2] = f[j+1, i, 2]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 6] = f[j+1, i+1, 6]

                f_new[j, i, 1] = f[j, i, 3]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 5] = f[j+1, i-1, 5]
                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 8] = f[j-1, i-1, 8]
            elif i == n_x - 1: # Top right corner
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 2] = f[j+1, i, 2]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 6] = f[j+1, i+1, 6]

                f_new[j, i, 1] = f[j, i, 3]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 5] = f[j+1, i-1, 5]
                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 8] = f[j-1, i-1, 8]
            else: # interior top surface
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 1] = f[j, i-1, 1]
                f_new[j, i, 2] = f[j+1, i, 2]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 5] = f[j+1, i-1, 5]
                f_new[j, i, 6] = f[j+1, i+1, 6]

                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 8] = f[j-1, i-1, 8]
        elif j == n_y - 1: # Bottom surface

            if i == 0: # Bottom left corner
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 1] = f[j, i-1, 1]
                f_new[j, i, 2] = f[j+1, i, 2]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 5] = f[j+1, i-1, 5]
                f_new[j, i, 6] = f[j+1, i+1, 6]
                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 8] = f[j-1, i-1, 8]

            if i == n_x - 1: # Bottom right corner
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 1] = f[j, i-1, 1]
                f_new[j, i, 2] = f[j+1, i, 2]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 5] = f[j+1, i-1, 5]
                f_new[j, i, 6] = f[j+1, i+1, 6]
                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 8] = f[j-1, i-1, 8]

            else: # interior bottom surface
                f_new[j, i, 0] = f[j, i, 0] 
                f_new[j, i, 1] = f[j, i-1, 1]
                f_new[j, i, 3] = f[j, i+1, 3]
                f_new[j, i, 4] = f[j-1, i, 4]
                f_new[j, i, 7] = f[j-1, i+1, 7]
                f_new[j, i, 8] = f[j-1, i-1, 8]

                f_new[j, i, 2] = f_new[j, i, 4]
                f_new[j, i, 5] = f_new[j, i, 7]
                f_new[j, i, 6] = f_new[j+1, i+1, 6]
        elif i == 0: # Left surface
            f_new[j, i, 0] = f[j, i, 0] 
            f_new[j, i, 2] = f[j+1, i, 2]
            f_new[j, i, 3] = f[j, i+1, 3]
            f_new[j, i, 4] = f[j-1, i, 4]
            f_new[j, i, 6] = f[j+1, i+1, 6]
            f_new[j, i, 7] = f[j-1, i+1, 7]

            # shorthand variables            
            f_0 = f_new[j, i, 0]
            f_2 = f_new[j, i, 2]
            f_3 = f_new[j, i, 3]
            f_4 = f_new[j, i, 4]
            f_6 = f_new[j, i, 6]
            f_7 = f_new[j, i, 7]

            # solve for unknown f
            U_in = 1 - (f_0+f_2+f_4+2*(f_3+f_6+f_7))/rho_in
            f_new[j, i, 1] = f[j, i, 3] + (2/3)*rho_in*U_in
            f_new[j, i, 5] = f_7 - 0.5*(f_2 - f_4) + (1/6)*rho_in*U_in
            f_new[j, i, 8] = f_6 + 0.5*(f_2 - f_4) + (1/6)*rho_in*U_in

        elif i == n_x - 1: # Right surface
            f_new[j, i, 0] = f[j, i, 0] 
            f_new[j, i, 1] = f[j, i-1, 1]
            f_new[j, i, 2] = f[j+1, i, 2]
            f_new[j, i, 4] = f[j-1, i, 4]
            f_new[j, i, 5] = f[j+1, i-1, 5]
            f_new[j, i, 8] = f[j-1, i-1, 8]

            f_new[j, i, 3] = f[j, i-1, 3]
            f_new[j, i, 6] = f[j, i-1, 6]
            f_new[j, i, 7] = f[j, i-1, 7]

        else: # Interior nodes
            f_new[j, i, 0] = f[j, i, 0] 
            f_new[j, i, 1] = f[j, i-1, 1]
            f_new[j, i, 2] = f[j+1, i, 2]
            f_new[j, i, 3] = f[j, i+1, 3]
            f_new[j, i, 4] = f[j-1, i, 4]
            f_new[j, i, 5] = f[j+1, i-1, 5]
            f_new[j, i, 6] = f[j+1, i+1, 6]
            f_new[j, i, 7] = f[j-1, i+1, 7]
            f_new[j, i, 8] = f[j-1, i-1, 8]
#Collisiop
for j in range(n_y):
    for i in range(n_x):
        # Compute moments
        rho[j, i] = np.sum(f[j, i])
        u_x[j, i] = np.sum(f[j, i]*ksi[:, 0])/rho[j, i]
        u_y[j, i] = np.sum(f[j, i]*ksi[:, 1])/rho[j, i]

        # 
        u = np.array([u_x[j, i], u_y[j, i]])
        f_eq[j, i] = fe_eq_formula(k, rho, u)

        f[j, i] = rho_in*weights
        rho[j, i] = np.sum(f[j, i]) + 1
#Update

#Visualization