import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

# open file given as argument
filename = sys.argv[1]
filenames = []

# create numpy tridimensional array to store velocity moduluses
all_U = []

with open(filename, 'r') as f:
    # read two integers: width and height of lattice
    line = f.readline()
    width, height = [int(x) for x in line.split()]
    length = width * height
    
    
    i = 0
    # while line is not empty
    while f.readline():
        # read length floats: densities
        densities = [float(x) for x in f.readline().split()]
        # read length Ux floats: x-velocities
        Ux = [float(x) for x in f.readline().split()]
        # read length Uy floats: y-velocities
        Uy = [float(x) for x in f.readline().split()]
        # compute velocity moduluses
        U = [np.sqrt(Ux[i]**2 + Uy[i]**2) for i in range(length)]
        
        # append velocity moduluses to array
        all_U.append(U)
        print(f'Frame {i}')
        i += 1

# Reshape the 1D list into a 2D array
all_U = np.array(all_U).reshape(-1, width, height)

# Function to update the plot for each frame
def update(frame):
    plt.clf()  # Clear the previous frame
    plt.imshow(all_U[frame], origin='upper', cmap='RdBu_r', vmin=0, vmax=0.2, interpolation='spline16')  # Adjust origin to top left and set color scale
    plt.title(f'Frame {frame}')
    plt.colorbar()
    
# Create the animation
animation = FuncAnimation(plt.figure(), update, frames=len(all_U), interval=100, repeat=False)

# Save the animation as a GIF
animation.save('movie.gif', writer='pillow')