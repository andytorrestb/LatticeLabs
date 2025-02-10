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

# Initialize computational domain
f_eq = np.zeros((n_y, n_x, 9))
f = np.zeros((n_y, n_x, 9))

for j in range(n_y):
    for i in range(n_x):
        for k in range(9):
            f_eq[j, i, k] = fe_eq_formula(k, rho[j, i], u_x[j, i], u_y[j, i])

            f = rho_in*weights[k]

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
            if i == n_x - 1: # Top right corner
            else: # interior top surfacs
        elif j == n_y - 1: # Bottom surface
            if i == 0: # Bottom left corner
            if i == n_x - 1: # Bottom right corner
            else: # interior bottom surface
        elif i == 0: # Left surface

        elif i == n_x - 1: # Right surface

        else: # Interior nodes
            f_new(j, i, 0) = f(j, i, 0) 
            f_new(j, i, 1) = f(j, i-1, 1)
            f_new(j, i, 2) = f(j+1, i, 2)
            f_new(j, i, 3) = f(j, i+1, 3)
            f_new(j, i, 4) = f(j-1, i, 4)
            f_new(j, i, 5) = f(j+1, i-1, 5)
            f_new(j, i, 6) = f(j+1, i+1, 6)
            f_new(j, i, 7) = f(j-1, i+1, 7)
            f_new(j, i, 8) = f(j-1, i-1, 8)
#Collisiop
#Update

#Visualization