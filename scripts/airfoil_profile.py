import numpy as np
import matplotlib.pyplot as plt
import os

def NACA_airfoil_profile(chord, height, thickness, camber, num_points=100):
    theta = np.linspace(0, 2 * np.pi, num_points)
    
    # NACA airfoil equation
    x = 0.5 * (1 - np.cos(theta)) * chord
    yt = 5 * thickness * (0.2969 * np.sqrt(x/chord) -
                         0.126 * (x/chord) - 0.3516 * (x/chord)**2 +
                         0.2843 * (x/chord)**3 - 0.1015 * (x/chord)**4)

    yc = camber * (1 - np.cos(theta))
    
    # Compute the y-coordinate of the airfoil profile above and below the camber line
    y_upper = yc + yt
    y_lower = yc - yt
    
    return x, y_upper, y_lower

def generate_grid(chord, height, thickness, camber, num_points=100, grid_size=200):
    grid = np.zeros((grid_size, grid_size), dtype=int)  # Initialize with 0
    
    x, y_upper, y_lower = NACA_airfoil_profile(chord, height, thickness, camber, num_points)
    
    # Normalize airfoil coordinates between -0.2 and 0.2 to center the wing in the grid
    x_normalized = (x - chord/2) / chord * 0.4
    y_upper_normalized = (y_upper - height/2) / height * 0.4
    y_lower_normalized = (y_lower - height/2) / height * 0.4
    
    # Find coordinates in the grid
    i_upper = np.round((y_upper_normalized + 0.5) * (grid_size - 1)).astype(int)
    i_lower = np.round((y_lower_normalized + 0.5) * (grid_size - 1)).astype(int)
    j = np.round((x_normalized + 0.5) * (grid_size - 1)).astype(int)
    
    # Set grid values above the airfoil profile line to 1
    grid[j, i_upper] = 1
    
    # Set grid values below the airfoil profile line to 1
    grid[j, i_lower] = 1
    
    return np.flipud(grid)  # Reverse the order of columns

def save_to_file_and_plot(grid, file_name):
    plt.imshow(grid, cmap='gray', origin='lower')
    plt.axis('off')
    plt.savefig(file_name)
    plt.close()

    with open(file_name.replace('.png', '.txt'), 'w') as file:
        for column in np.flipud(grid.T):  # Reverse the order of columns
            file.write(' '.join(map(lambda x: '1' if x == 0 else '0', column)) + '\n')

def replace_1_with_0(line):
    for i in range(len(line)):
        if line[i] == '1':
            if '0' in line[:i] and '0' in line[i+1:]:
                line = line[:i] + '0' + line[i+1:]
    return line

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    modified_lines = [replace_1_with_0(line) for line in lines]

    with open(output_file, 'w') as outfile:
        outfile.writelines(modified_lines)

def plot_from_file(file_name):
    modified_grid = np.loadtxt(file_name, dtype=int)
    plt.imshow(modified_grid, cmap='gray', origin='lower')
    plt.axis('off')

    # Save the image to a PNG file
    image_file_name = file_name.replace('.txt', '_plot.png')
    plt.savefig(image_file_name)

    print(f"Image saved to '{image_file_name}'.")

def read_file2(file_path):
    with open(file_path, 'r') as file:
        return [list(line.strip()) for line in file]

def rotate_clockwise(matrix):
    # Rotate the matrix clockwise
    return [list(row) for row in zip(*reversed(matrix))]

def write_file2(file_path, matrix):
    with open(file_path, 'w') as file:
        for row in matrix:
            file.write(''.join(row) + '\n')

def main2():
    # Enter the path to your text file
    input_file_path = 'airfoil_profile_new.txt'
    
    # Read the matrix from the text file input
    original_matrix = read_file2(input_file_path)

    # Rotate the image clockwise
    rotated_matrix = rotate_clockwise(original_matrix)

    # Enter the path for the new output file
    output_file_path = 'airfoil_profile_new_rotated.txt'

    # Write the rotated image to the new file
    write_file2(output_file_path, rotated_matrix)

def read_file3(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def write_file3(file_path, data):
    inverted_data = [''.join(['1' if c == '0' else '0' for c in row]) for row in data]
    with open(file_path, 'w') as file:
        file.write('\n'.join(inverted_data))

def display_plot3(data, output_path):
    matrix_np = np.array([[0 if c == '0' else 1 for c in row] for row in data])
    plt.imshow(matrix_np, cmap='gray', interpolation='none')
    plt.savefig(output_path, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.show()

def main3():
    # Enter the path to your text file
    input_file_path = 'airfoil_profile_new_rotated.txt'
    
    # Read the data from the text file, removing empty lines
    data = read_file3(input_file_path)

    # Enter the path for the new output file
    output_file_path = 'output_plot.png'

    # Write the new modified file by inverting zeros and ones
    new_modified_file_path = 'new_modified_file.txt'
    write_file3(new_modified_file_path, data)

    # Display the plot with black 0s and white 1s
    display_plot3(data, output_file_path)


if __name__ == "__main__":
    chord = 0.7        # Wing chord length
    height = 0.06     # Wing height
    thickness = 0.01  # Wing thickness
    camber = 0.01      # Wing camber
    
    num_points = 500   # Number of points on the wing chord
    grid_size = 300    # Grid size
    
    # Generate and save the original grid
    grid = generate_grid(chord, height, thickness, camber, num_points, grid_size)
    file_name = 'airfoil_profile.png'
    save_to_file_and_plot(grid, file_name)
    print(f"Result saved to '{file_name}'.")

    # Modify the text file
    input_file_name = "airfoil_profile.txt"
    output_file_name = "airfoil_profile_new.txt"
    process_file(input_file_name, output_file_name)

    # Plot the image from the modified file
    plot_from_file(output_file_name)
    main2()
    main3()

    generated_files = os.listdir()
for file_name in generated_files:
    if file_name.endswith('.txt') or file_name.endswith('.png'):
        if file_name not in ['output_plot.png', 'new_modified_file.txt']:
            os.remove(file_name)
