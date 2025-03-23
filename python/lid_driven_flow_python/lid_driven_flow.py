import numpy as np

# Parameters
N_x = 100
N_y = 50
dx, dy = 1, 1
Ksi = np.array([[0, 1, 0, -1, 0, 1, -1, -1, 1],
                [0, 0, 1, 0, -1, 1, 1, -1, -1]])
w = np.array([4/9] + [1/9]*4 + [1/36]*4)
c_s = 1/np.sqrt(3)
Tau = 1

# Initialization
rho_in = 1
rho = np.ones((N_y, N_x)) * rho_in
u = np.zeros((N_y, N_x))
v = np.zeros((N_y, N_x))
f_eq = np.zeros((N_y, N_x, 9))
f = np.zeros((N_y, N_x, 9))

for j in range(N_y):
    for i in range(N_x):
        for k in range(9):
            vel = np.dot(Ksi[:, k], [u[j, i], v[j, i]])
            f_eq[j, i, k] = w[k] * rho[j, i] * (1 + vel / c_s**2 + vel**2 / (2*c_s**4) - (u[j, i]**2 + v[j, i]**2) / (2*c_s**2))
            f[j, i, k] = rho_in * w[k]

f_new = np.zeros((N_y, N_x, 9))

T = 20000
for t in range(T):
    print('Time step:', t)
    # Streaming
    for j in range(N_y):
        for i in range(N_x):
            if j == 0:  # Top surface
                if i == 0:  # Top-Left corner
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 2] = f[j + 1, i, 2]
                    f_new[j, i, 3] = f[j, i + 1, 3]
                    f_new[j, i, 6] = f[j + 1, i + 1, 6]

                    f_new[j, i, 1] = f_new[j, i, 3]
                    f_new[j, i, 4] = f_new[j, i, 2]
                    f_new[j, i, 5] = rho_in / 2 - f_new[j, i, 0] / 2 - f_new[j, i, 2] - f_new[j, i, 3] - f_new[j, i, 6]
                    f_new[j, i, 7] = f_new[j, i, 5]
                    f_new[j, i, 8] = f_new[j, i, 6]

                elif i == N_x - 1:  # Top-Right corner
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 1] = f[j, i - 1, 1]
                    f_new[j, i, 2] = f[j + 1, i, 2]
                    f_new[j, i, 5] = f[j + 1, i - 1, 5]

                    f_new[j, i, 3] = f_new[j, i, 1]
                    f_new[j, i, 4] = f_new[j, i, 2]
                    f_new[j, i, 6] = f_new[j, i - 1, 6]
                    f_new[j, i, 7] = f_new[j, i, 5]
                    f_new[j, i, 8] = f_new[j, i, 6]

                else:  # All other nodes on the top surface
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 1] = f[j, i - 1, 1]
                    f_new[j, i, 2] = f[j + 1, i, 2]
                    f_new[j, i, 3] = f[j, i + 1, 3]
                    f_new[j, i, 5] = f[j + 1, i - 1, 5]
                    f_new[j, i, 6] = f[j + 1, i + 1, 6]

                    f_new[j, i, 4] = f_new[j, i, 2]
                    f_new[j, i, 7] = (f_new[j, i, 1] - f_new[j, i, 3]) / 2 + f_new[j, i, 5]
                    f_new[j, i, 8] = (f_new[j, i, 3] - f_new[j, i, 1]) / 2 + f_new[j, i, 6]

            elif j == N_y - 1:  # Bottom surface
                if i == 0:  # Bottom-Left corner
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 3] = f[j, i + 1, 3]
                    f_new[j, i, 4] = f[j - 1, i, 4]
                    f_new[j, i, 7] = f[j - 1, i + 1, 7]

                    f_new[j, i, 1] = f_new[j, i, 3]
                    f_new[j, i, 2] = f_new[j, i, 4]
                    f_new[j, i, 5] = f_new[j, i, 7]
                    f_new[j, i, 6] = rho_in / 2 - f_new[j, i, 0] / 2 - f_new[j, i, 3] - f_new[j, i, 4] - f_new[j, i, 7]
                    f_new[j, i, 8] = f_new[j, i, 6]

                elif i == N_x - 1:  # Bottom-Right corner
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 1] = f[j, i - 1, 1]
                    f_new[j, i, 4] = f[j - 1, i, 4]
                    f_new[j, i, 8] = f[j - 1, i - 1, 8]

                    f_new[j, i, 2] = f_new[j, i, 4]
                    f_new[j, i, 3] = f_new[j, i, 1]
                    f_new[j, i, 5] = f_new[j, i - 1, 5]
                    f_new[j, i, 6] = f_new[j, i, 8]
                    f_new[j, i, 7] = f_new[j, i, 5]

                else:  # All other nodes on the bottom surface
                    f_new[j, i, 0] = f[j, i, 0]
                    f_new[j, i, 1] = f[j, i - 1, 1]
                    f_new[j, i, 3] = f[j, i + 1, 3]
                    f_new[j, i, 4] = f[j - 1, i, 4]
                    f_new[j, i, 7] = f[j - 1, i + 1, 7]
                    f_new[j, i, 8] = f[j - 1, i - 1, 8]

                    f_new[j, i, 2] = f_new[j, i, 4]
                    f_new[j, i, 5] = (f_new[j, i, 3] - f_new[j, i, 1]) / 2 + f_new[j, i, 7]
                    f_new[j, i, 6] = (f_new[j, i, 1] - f_new[j, i, 3]) / 2 + f_new[j, i, 8]

            elif i == 0:  # Left surface
                f_new[j, i, 0] = f[j, i, 0]
                f_new[j, i, 2] = f[j + 1, i, 2]
                f_new[j, i, 3] = f[j, i + 1, 3]
                f_new[j, i, 4] = f[j - 1, i, 4]
                f_new[j, i, 6] = f[j + 1, i + 1, 6]
                f_new[j, i, 7] = f[j - 1, i + 1, 7]

                U_in = 1 - (f_new[j, i, 0] + f_new[j, i, 2] + f_new[j, i, 4] + 2 * (f_new[j, i, 3] + f_new[j, i, 6] + f_new[j, i, 7])) / rho_in
                f_new[j, i, 1] = f_new[j, i, 3] + U_in * rho_in * 2 / 3
                f_new[j, i, 5] = f_new[j, i, 7] + (f_new[j, i, 4] - f_new[j, i, 2]) / 2 + U_in * rho_in / 6
                f_new[j, i, 8] = f_new[j, i, 6] - (f_new[j, i, 4] - f_new[j, i, 2]) / 2 + U_in * rho_in / 6

            elif i == N_x - 1:  # Right surface
                f_new[j, i, 0] = f[j, i, 0]
                f_new[j, i, 1] = f[j, i - 1, 1]
                f_new[j, i, 2] = f[j + 1, i, 2]
                f_new[j, i, 4] = f[j - 1, i, 4]
                f_new[j, i, 5] = f[j + 1, i - 1, 5]
                f_new[j, i, 8] = f[j - 1, i - 1, 8]

                f_new[j, i, 3] = f_new[j, i - 1, 3]
                f_new[j, i, 6] = f_new[j, i - 1, 6]
                f_new[j, i, 7] = f_new[j, i - 1, 7]

            else:  # All interior nodes
                f_new[j, i, 0] = f[j, i, 0]
                f_new[j, i, 1] = f[j, i - 1, 1]
                f_new[j, i, 2] = f[j + 1, i, 2]
                f_new[j, i, 3] = f[j, i + 1, 3]
                f_new[j, i, 4] = f[j - 1, i, 4]
                f_new[j, i, 5] = f[j + 1, i - 1, 5]
                f_new[j, i, 6] = f[j + 1, i + 1, 6]
                f_new[j, i, 7] = f[j - 1, i + 1, 7]
                f_new[j, i, 8] = f[j, i - 1, 8]



    # Collision and update moments
    for j in range(N_y):
        for i in range(N_x):
            rho[j, i] = np.sum(f_new[j, i, :])
            u[j, i] = (f_new[j, i, 1] + f_new[j, i, 5] + f_new[j, i, 8] - f_new[j, i, 3] - f_new[j, i, 6] - f_new[j, i, 7]) / rho[j, i]
            v[j, i] = (f_new[j, i, 2] + f_new[j, i, 5] + f_new[j, i, 6] - f_new[j, i, 4] - f_new[j, i, 7] - f_new[j, i, 8]) / rho[j, i]

    for j in range(N_y):
        for i in range(N_x):
            vel = np.clip(np.dot(Ksi.T, [u[j, i], v[j, i]]), -0.2*c_s, 0.2*c_s)
            f_eq[j, i, :] = rho[j, i] * w * (1 + vel / c_s**2 + vel**2 / (2*c_s**4) - (u[j, i]**2 + v[j, i]**2) / (2*c_s**2))

    f = f_new - (1/Tau) * (f_new - f_eq)




import numpy as np
import matplotlib.pyplot as plt

# Benchmark calculation
y_benchmark = np.arange(0, 1.01, 0.01)
u_benchmark = -4 * (y_benchmark**2 - y_benchmark)

# Plot benchmark result
plt.figure()
plt.plot(y_benchmark, u_benchmark, "r", label="Benchmark")

# Simulation result (assuming u and N_x are defined)
u_sim = np.zeros(N_y-1)  # Placeholder for actual u_sim computation
print(u.shape)
print(N_y)
for j in range(N_y-1):
    print(j)
    u_sim[j] = u[j, N_x - 1]  # Ensure u and N_x are defined

plt.plot(np.linspace(0, 1, N_y-1), u_sim / np.max(u_sim), "b", label="Simulation")
plt.legend()
plt.savefig('velocity_profile.png')

# Quiver plot (assuming u and v are defined)
plt.figure()
plt.quiver(np.flipud(u), np.flipud(v), scale=10)
plt.axis("equal")
plt.axis("tight")
plt.savefig('velocity.png')

# Contour plot (assuming Rho is defined)
plt.figure()
plt.contourf(rho, levels=30)
plt.axis("equal")
plt.axis("tight")
plt.colorbar()
plt.savefig('density.png')
