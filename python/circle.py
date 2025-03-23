import random
import numpy as np
import matplotlib.pyplot as plt

def generate_a_point():
    x = np.round(random.random(), 3)
    y = np.round(random.random(), 3)
    return x, y

def test_circle(x, y, R, center_x, center_y):
    return (x - center_x)**2 + (y - center_y)**2 <= R**2

def generate_an_inside_point(R, center_x, center_y):
    x = np.round(random.uniform(center_x - R, center_x + R), 3)
    y = np.round(random.uniform(center_y - R, center_y + R), 3)

    if test_circle(x, y, R, center_x, center_y):
        return x, y
    else:
        return generate_an_inside_point(R, center_x, center_y)


def generate_an_outside_point(R, center_x, center_y):

    x_bounds = (-5, 5)
    y_bounds = (-5, 5)

    x_min, y_min = x_bounds[0], y_bounds[0]
    x_max, y_max = x_bounds[1], y_bounds[1]

    min_val = min(x_min, y_min)
    max_val = max(x_max, y_max)


    x = np.round(random.uniform(x_min, x_max), 3)
    y = np.round(random.uniform(y_min, y_max), 3)

    if test_circle(x, y, R, center_x, center_y):
        return generate_an_outside_point(R, center_x, center_y)
    else:
        return x, y


def test_inside_points():
    R = 1
    center_x = 0.0
    center_y = 0.0

    # Plot the circle
    circle = plt.Circle((center_x, center_y), R, color='blue', fill=False)
    fig, ax = plt.subplots()
    ax.add_artist(circle)

    for i in range(100):
        x, y = generate_an_inside_point(R, center_x, center_y)
        plt.plot(x, y, 'ro' if test_circle(x, y, R, center_x, center_y) else 'go')

    ax.set_aspect('equal', 'box')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)

    # Show the plot
    plt.savefig('inside_circle.png')

def test_outside_points():
    R = 1
    center_x = 0.0
    center_y = 0.0

    # Plot the circle
    circle = plt.Circle((center_x, center_y), R, color='blue', fill=False)
    fig, ax = plt.subplots()
    ax.add_artist(circle)

    for i in range(100):
        x, y = generate_an_outside_point(R, center_x, center_y)
        plt.plot(x, y, 'ro' if test_circle(x, y, R, center_x, center_y) else 'go')

    ax.set_aspect('equal', 'box')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)

    # Show the plot
    plt.savefig('outside_circle.png')

test_inside_points()
test_outside_points()


def calculate_intersection():

    # Define circle
    center_x, center_y = 0.0, 0.0
    R = 1.0

    x_in, y_in = generate_an_inside_point(R, center_x, center_y)
    x_out, y_out = generate_an_outside_point(R, center_x, center_y)

    m = (y_out - y_in) / (x_out - x_in)

    x, y = x_in, y_in
    R_point = np.sqrt(x**2+y**2)

    dx = 0.01

    while R_point < R:
        x = x + dx
        y = m*x + b

    return x,y


    while 
# # Generate a point and test if it's inside the circle
# x, y = generate_a_point()
# R = 1
# center_x = 0.0
# center_y = 0.0
# is_inside = test_circle(x, y, R, center_x, center_y)

# print(f'Point ({x}, {y}) is inside the circle: {is_inside}')
# # Plot the circle
# circle = plt.Circle((center_x, center_y), R, color='blue', fill=False)
# fig, ax = plt.subplots()
# ax.add_artist(circle)

# # Plot the point
# ax.plot(x, y, 'ro' if is_inside else 'go')  # Red if inside, green if outside

# # Set the aspect of the plot to be equal
# ax.set_aspect('equal', 'box')
# ax.set_xlim(-1.5, 1.5)
# ax.set_ylim(-1.5, 1.5)

# # Show the plot
# plt.savefig('circle.png')
