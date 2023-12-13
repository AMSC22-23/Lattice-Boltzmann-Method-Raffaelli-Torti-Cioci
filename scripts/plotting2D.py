import matplotlib.pyplot as plt
import sys
import numpy as np
import os
import imageio

# open file given as argument
filename = sys.argv[1]
filenames = []

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
        
        # create animation frame with velocity moduluses only
        plt.figure(figsize=(width/10, height/10))
        plt.title('Velocity moduluses')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(0, width)
        plt.ylim(0, height)
        plt.pcolor(np.reshape(U, (height, width)))
        plt.colorbar()
        plt.savefig('frame' + str(i) + '.png')
        plt.close()
        i += 1
        filenames.append('frame' + str(i) + '.png')
        
# remove last filename
filenames.pop()
# create gif from png with python code
with imageio.get_writer('movie.gif', mode='I') as writer:
    for f in filenames:
        image = imageio.imread(f)
        writer.append_data(image)

# remove frames    
os.system('rm frame*.png')