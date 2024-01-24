import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath

def naca4(number, n=100, tilt_angle=0):
    """ Create the NACA 4-digit airfoil profile """
    m = int(number[0]) / 100.0
    p = int(number[1]) / 10.0
    t = int(number[2:]) / 100.0

    x = np.linspace(0, 1, n)
    yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

    if p != 0:
        yc = np.where(x <= p, m / p**2 * (2*p*x - x**2), m / (1-p)**2 * ((1 - 2*p) + 2*p*x - x**2))
        dyc_dx = np.where(x <= p, 2*m / p**2 * (p - x), 2*m / (1-p)**2 * (p - x))
        theta = np.arctan(dyc_dx)
    else:
        yc = np.zeros_like(x)
        theta = np.zeros_like(x)

    xu = x - yt*np.sin(theta)
    xl = x + yt*np.sin(theta)
    yu = yc + yt*np.cos(theta)
    yl = yc - yt*np.cos(theta)

    # Rotate the airfoil
    tilt_angle_rad = np.deg2rad(tilt_angle)
    cos_angle = np.cos(tilt_angle_rad)
    sin_angle = np.sin(tilt_angle_rad)
    xu, yu = cos_angle * xu - sin_angle * yu, sin_angle * xu + cos_angle * yu
    xl, yl = cos_angle * xl - sin_angle * yl, sin_angle * xl + cos_angle * yl

    return xu, yu, xl, yl
            
def discretize_airfoil(xu, yu, xl, yl, grid_width, grid_height, padding_left, padding_right, padding_vertical):
    """ Discretize the airfoil onto a grid of 0s and 1s """
    # Create the airfoil path
    airfoil_path = mpath.Path(np.vstack((np.hstack((xu, xl[::-1])), np.hstack((yu, yl[::-1])))).T)

    # Create the grid
    x = np.linspace(np.min(xu) - padding_left, np.max(xu) + padding_right, grid_width)
    y = np.linspace(np.min([yu, yl]) - padding_vertical, np.max([yu, yl]) + padding_vertical, grid_height)
    X, Y = np.meshgrid(x, y)

    # Check if each point in the grid is inside the airfoil
    grid = airfoil_path.contains_points(np.vstack((X.flatten(), Y.flatten())).T).reshape(grid_height, grid_width)

    return grid

def write_airfoil_to_grid(xu, yu, xl, yl, grid_width, grid_height, padding_left, padding_right, padding_vertical, filename):
    """ Write the grid coordinates of the airfoil to a text file """
    # Discretize the airfoil
    grid = discretize_airfoil(xu, yu, xl, yl, grid_width, grid_height, padding_left, padding_right, padding_vertical)

    # Create the grid coordinates
    x = np.linspace(np.min(xu) - padding_left, np.max(xu) + padding_right, grid_width)
    y = np.linspace(np.min([yu, yl]) - padding_vertical, np.max([yu, yl]) + padding_vertical, grid_height)
    X, Y = np.meshgrid(x, y)

    # Write the grid coordinates of the airfoil to the file
    with open(filename, 'w') as f:
        for i in range(grid_height):
            for j in range(grid_width):
                if grid[i, j]:
                    f.write(f'{j} {grid_height - i - 1}\n')  # Flip the y-coordinate to match the image orientation


def plot_airfoil_from_file(filename, grid_width, grid_height):
    """ Read the grid coordinates of the airfoil from a text file and plot them """
    # Create an empty grid
    grid = np.zeros((grid_height, grid_width))

    # Read the grid coordinates from the file
    with open(filename, 'r') as f:
        for line in f:
            x, y = map(int, line.split())
            grid[grid_height - y - 1, x] = 1  # Flip the y-coordinate to match the image orientation

    # Plot the grid
    plt.imshow(grid, cmap='gray_r', origin='lower')
    plt.savefig('airfoil_from_file.png')

# change parameters here
xu, yu, xl, yl = naca4('2412', tilt_angle=-15)
write_airfoil_to_grid(xu, yu, xl, yl, grid_width=400, grid_height=150, padding_left=0.1, padding_right=0.6, padding_vertical=0.3, filename='airfoil_coordinates.txt')
plot_airfoil_from_file('airfoil_coordinates.txt', grid_width=400, grid_height=150)
