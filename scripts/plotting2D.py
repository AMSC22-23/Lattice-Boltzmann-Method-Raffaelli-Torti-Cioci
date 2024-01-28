import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import time

# open file given as argument
filename = sys.argv[1]

# create numpy tridimensional array to store velocity moduluses
all_U = []
all_Ux = []
all_Uy = []
all_Steps = []

with open(filename, 'r') as f:
    # read two integers: width and height of lattice
    line = f.readline()
    width, height = [int(x) for x in line.split()]
    length = width * height
    
    i = 0
    # while line is not empty
    line = f.readline()
    while line:
        step = line
        # read length Ux floats: x-velocities
        line = f.readline()  # read next line
        Ux = [float(x) for x in line.split()]
        # read length Uy floats: y-velocities
        line = f.readline()  # read next line
        Uy = [float(x) for x in line.split()]
        # compute velocity moduluses
        U = [np.sqrt(Ux[i]**2 + Uy[i]**2) for i in range(length)]
        
        all_U.append(U)
        all_Ux.append(Ux)
        all_Uy.append(Uy)
        all_Steps.append(step)
        print(f'Frame {i}')
        i += 1
        line = f.readline()

# Reshape the 1D list into a 2D array
all_U = np.array(all_U).reshape(-1, height, width)
all_Ux = np.array(all_Ux).reshape(-1, height, width)
all_Uy = np.array(all_Uy).reshape(-1, height, width)

# Function to update the plot for each frame
def update(frame):
    plt.clf()  # Clear the previous frame
    
    dfydx = all_Ux[frame, 2:, 1:-1] - all_Ux[frame, :-2, 1:-1]
    dfxdy = all_Uy[frame, 1:-1, 2:] - all_Uy[frame, 1:-1, :-2]
    curl = dfydx - dfxdy
    
    plt.imshow(all_U[frame], origin='upper', cmap='RdBu_r', vmin=0, vmax=0.3, interpolation='spline16')
    plt.colorbar()
    plt.title(f'Step {all_Steps[frame]}')
    
    '''
    # Calculate the step size for x and y directions
    step_x = all_Ux.shape[2] // 6
    step_y = all_Ux.shape[1] // 6

    # Create the meshgrid
    X, Y = np.meshgrid(np.arange(step_x//2, all_Ux.shape[2], step_x), np.arange(step_y//2, all_Ux.shape[1], step_y))

    # Plot the quiver plot
    plt.quiver(X, Y, all_Ux[frame, Y, X], all_Uy[frame, Y, X], color='white', angles='xy', scale=2.5)
    '''

# Create the animation
animation = FuncAnimation(plt.figure(), update, frames=len(all_U), interval=100, repeat=False)

# Save the animation as video with system current time in the name
timestr = time.strftime("%Y%m%d-%H%M%S")
animation.save(f'outputs/movie-{timestr}.mp4', fps=10, writer='ffmpeg')