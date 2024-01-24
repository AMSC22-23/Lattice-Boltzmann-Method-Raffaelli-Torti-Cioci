import numpy as np
import matplotlib.pyplot as plt
import os

# Define the folder name
folder_name = 'src'

# Get the current working directory
current_directory = os.getcwd()

# Construct the full path to the data file
file_path = os.path.join(current_directory, folder_name, 'dragLift.txt')

# Read data from the dragLift.txt file
data = np.loadtxt(file_path, dtype=float)

# Extract columns for time instant, drag, and lift
time_instant = data[:, 0]
drag = data[:, 1]
lift = data[:, 2]

# Create the first plot for drag
plt.figure(figsize=(10, 5))
plt.subplot(2, 1, 1)
plt.plot(time_instant, drag, label='Drag', color='red')
plt.title('Drag Plot')
plt.xlabel('Time Instant')
plt.ylabel('Drag')
plt.legend()

# Create the second plot for lift
plt.subplot(2, 1, 2)
plt.plot(time_instant, lift, label='Lift', color='blue')
plt.title('Lift Plot')
plt.xlabel('Time Instant')
plt.ylabel('Lift')
plt.legend()

# Show the plots
plt.tight_layout()
plt.show()
