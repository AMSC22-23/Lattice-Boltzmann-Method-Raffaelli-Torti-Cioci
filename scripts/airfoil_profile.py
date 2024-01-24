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

    with open(file_name.replace('_AF.png', '_AF.txt'), 'w') as file:
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
    image_file_name = file_name.replace('_AF.txt', '_plot_AF.png')
    plt.savefig(image_file_name)

    print(f"Image saved to '{image_file_name}'.")

    

def main1():
    chord = 0.8        # Wing chord length
    height = 0.06     # Wing height
    thickness = 0.01  # Wing thickness
    camber = 0.03      # Wing camber
    
    num_points = 500   # Number of points on the wing chord
    grid_size = 300    # Grid size
    
    # Generate and save the original grid
    grid = generate_grid(chord, height, thickness, camber, num_points, grid_size)
    file_name = 'airfoil_profile_AF.png'
    save_to_file_and_plot(grid, file_name)
    print(f"Result saved to '{file_name}'.")

    # Modify the text file
    input_file_name = "airfoil_profile_AF.txt"
    output_file_name = "airfoil_profile_new_AF.txt"
    process_file(input_file_name, output_file_name)

    # Plot the image from the modified file
    plot_from_file(output_file_name)





    

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
    input_file_path = 'airfoil_profile_new_AF.txt'
    
    # Read the matrix from the text file input
    original_matrix = read_file2(input_file_path)

    # Rotate the image clockwise
    rotated_matrix_t = rotate_clockwise(original_matrix)

    rotated_matrix = rotate_clockwise(rotated_matrix_t)

    # Enter the path for the new output file
    output_file_path = 'airfoil_profile_new_rotated_AF.txt'

    # Write the rotated image to the new file
    write_file2(output_file_path, rotated_matrix)





def read_file3(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def write_file3(file_path, content):
    with open(file_path, 'w') as file:
        for line in content:
            modified_line = line.replace('0', 'x').replace('1', '0').replace('x', '1')
            file.write(modified_line + '\n')


            

def plot_binary_file(file_path, output_image_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]

    # Height and width of the image in pixels
    height = len(lines)
    width = len(lines[0])

    # Create a 2D array of zeros
    image_array = [[0] * width for _ in range(height)]

    # Fill the array with file values, skip non-'0' and non-'1' characters
    for i, line in enumerate(lines):
        for j, char in enumerate(line):
            if char in ['0', '1']:
                image_array[i][j] = int(char)
            elif char != ' ':
                print(f"Warning: Skipping invalid character '{char}' at position ({i}, {j})")

    # Create the plot
    fig, ax = plt.subplots()

    # Draw white pixels (0) and black pixels (1)
    ax.imshow(image_array, cmap='gray', interpolation='none')

    # Set the aspect ratio of the image to be equal
    ax.set_aspect('equal')

    # Remove axis labels
    ax.set_xticks([])
    ax.set_yticks([])

    # Save the image
    plt.savefig(output_image_path, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.show()



def extract_coordinates(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = [line.replace(' ', '') for line in file]

    coordinates = []

    for i, line in enumerate(lines):
        for j, char in enumerate(line):
            if char == '1':
                coordinates.append((i, j))

    with open(output_file, 'w') as output:
        for coord in coordinates:
            output.write(f"{coord[0]} {coord[1]}\n")


def transform_coordinates(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read coordinates from the first line
            first_line = infile.readline().strip().split()
            if len(first_line) != 2:
                raise ValueError("The file must contain pairs of x and y coordinates separated by a space.")
            
            # Calculate the transformation to apply to other coordinates
            x_offset, y_offset = float(first_line[0]), float(first_line[1])

            # Write the transformed first pair of coordinates (rounded to integers) to the new file
            outfile.write("0 0\n")

            # Read and transform the remaining pairs of coordinates
            for line in infile:
                coordinates = line.strip().split()
                if len(coordinates) != 2:
                    raise ValueError("The file must contain pairs of x and y coordinates separated by a space.")
                
                x, y = float(coordinates[0]), float(coordinates[1])
                x_transformed = round(x - x_offset)
                y_transformed = round(y - y_offset)

                # Write the transformed coordinates (rounded to integers) to the new file
                outfile.write(f"{x_transformed} {y_transformed}\n")

        print(f"Transformation completed. The result has been written to '{output_file}'.")
    
    except Exception as e:
        print(f"An error occurred: {e}")



def add_offset(input_file, output_file, row, column):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read and add the offset to each pair of coordinates
            for line in infile:
                coordinates = line.strip().split()
                if len(coordinates) != 2:
                    raise ValueError("The file must contain pairs of x and y coordinates separated by a space.")

                x, y = int(coordinates[0]) + column, int(coordinates[1]) + row

                # Write the pair of coordinates with the offset to the new file
                outfile.write(f"{x} {y}\n")

        print(f"Offset addition completed. The result has been written to '{output_file}'.")
    
    except Exception as e:
        print(f"An error occurred: {e}")


def rotate_coordinates(x, y, angle):
    """
    Ruota le coordinate (x, y) di un angolo specificato mantenendo le coordinate come interi positivi.
    """
    angle_rad = np.radians(angle)
    x_rotated = int(np.round(x * np.cos(angle_rad) - y * np.sin(angle_rad)))
    y_rotated = int(np.round(x * np.sin(angle_rad) + y * np.cos(angle_rad)))
    return x_rotated, y_rotated

def rotate_and_save_coordinates(input_file, output_file, angle):
    """
    Legge le coppie di coordinate dal file di input, ruota le coordinate e salva il risultato nel file di output.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            x, y = map(int, line.strip().split())
            x_rotated, y_rotated = rotate_coordinates(x, y, angle)
            outfile.write(f"{x_rotated} {y_rotated}\n")


def create_binary_file(input_file, output_file):
    # Inizializza un dizionario per tenere traccia delle coordinate esistenti
    coordinates = {}

    # Leggi il file di input e registra le coordinate nel dizionario
    with open(input_file, 'r') as f:
        for line in f:
            x, y = map(int, line.strip().split())
            coordinates[(x, y)] = 1

    # Scrivi nel file di output con zeri e uni in base alla presenza delle coordinate
    with open(output_file, 'w') as f:
        for x in range(400):  # Puoi regolare il range in base alle tue coordinate massime
            for y in range(400):
                if (x, y) in coordinates:
                    f.write('1 ')
                else:
                    f.write('0 ')
            f.write('\n')


def replace_zeros_between_ones(line):
    found_first_one = False
    modified_line = []

    for char in line:
        if char == '1':
            found_first_one = True
            modified_line.append(char)
        elif found_first_one and char == '0':
            modified_line.append('1')
        else:
            modified_line.append(char)

    return ''.join(modified_line)

def replace_zeros_between_ones(line):
    found_first_one = False
    modified_line = []

    for char in line:
        if char == '1':
            found_first_one = True
            modified_line.append(char)
        elif found_first_one and char == '0':
            modified_line.append('1')
        else:
            modified_line.append(char)

    return ''.join(modified_line)

def replace_zeros_between_ones(line):
    found_first_one = False
    modified_line = []

    for char in line:
        if char == '1':
            found_first_one = True
            modified_line.append(char)
        elif found_first_one and char == '0':
            # Smetti di sostituire '0' con '1' quando trovi l'ultimo '1'
            if '1' in modified_line:
                found_first_one = False
            modified_line.append('1')
        else:
            modified_line.append(char)

    return ''.join(modified_line)

def process_file_2(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Crea il file di output se non esiste
    if not os.path.exists(output_file):
        open(output_file, 'w').close()

    modified_lines = [replace_zeros_between_ones(line.strip()) for line in lines]

    with open(output_file, 'w') as outfile:
        for modified_line in modified_lines:
            outfile.write(modified_line + '\n')


def plot_and_save_image(file_path, output_path='output.png'):
    try:
        # Read the file of 1s and 0s
        with open(file_path, 'r') as file:
            # Read all lines and remove spaces
            lines = [line.replace(' ', '').strip() for line in file]

        # Remove any empty lines
        lines = [line for line in lines if line]

        # Transform the string into a list of lists of integers
        image_data = [[int(char) for char in row] for row in lines]

        # Create a NumPy array for the plot
        image_array = np.array(image_data)

        # Create the plot
        plt.imshow(image_array, cmap='gray', interpolation='nearest')

        # Set colors for 0 and 1
        plt.colorbar(ticks=[0, 1])

        # Set axis labels
        plt.xlabel('Column')
        plt.ylabel('Row')

        # Save the plot as a PNG file
        plt.savefig(output_path, format='png')

        # Display the plot
        plt.show()

        print(f"Image saved as '{output_path}'")

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unknown error: {e}")
















def main3():
    # Enter the path to your text file
    input_file_path = 'airfoil_profile_new_rotated_AF.txt'
    
    # Read the data from the text file, removing empty lines
    data = read_file3(input_file_path)

    # Enter the path for the new output file
    output_image_path = 'airfoil_plot_AF.png'

    # Write the new modified file by inverting zeros and ones
    new_modified_file_path = 'visual_file_AF.txt'
    write_file3(new_modified_file_path, data)

    # Display the plot with black 0s and white 1s
    # plot_from_file(new_modified_file_path)

    # Plot the binary file and save the image
    # plot_binary_file(new_modified_file_path, output_image_path)




    # get coordinates
    input_file = 'visual_file_AF.txt'  # Replace with the actual path of your input file
    output_file = 'airfoil_coordinates_AF.txt'  # Replace with the desired output file path

    extract_coordinates(input_file, output_file)


    input_file = 'airfoil_coordinates_AF.txt'  # Replace with the actual path of your input file
    output_file = 'airfoil_coordinates_no_offset_AF.txt'  # Replace with the desired output file path

    # Example usage:
    transform_coordinates(input_file, output_file)

    input_file = 'airfoil_coordinates_no_offset_AF.txt'
    output_file = 'airfoil_coordinates_no_offset_angle_AF.txt'
    angle = 65  # Angolo di rotazione in gradi

    rotate_and_save_coordinates(input_file, output_file, angle)


    input_file = 'airfoil_coordinates_no_offset_angle_AF.txt'  # Replace with the actual path of your input file
    output_file = 'airfoil_coordinates_offset_AF.txt'  # Replace with the desired output file path
    # Example usage:
    add_offset(input_file, output_file, 150, 150)


    input_file = "airfoil_coordinates_offset_AF.txt"
    output_file = "visual_not_processed_AF.txt"
    create_binary_file(input_file, output_file)

    input_file = "visual_not_processed_AF.txt"
    output_file = "visual_file_rotated_AF.txt"
    process_file_2(input_file, output_file)
    extract_coordinates(output_file, "airfoil_coordinates_offset_AF.txt")



    # Esempio di utilizzo con un percorso di output personalizzato
    file_path = 'visual_file_rotated_AF.txt'
    output_path = 'airfoil_plot_AF.png'
    plot_and_save_image(file_path, output_path)






if __name__ == "__main__":
    main1()
    main2()
    main3()


    generated_files = os.listdir()
    for file_name in generated_files:
        if file_name.endswith('_AF.txt') or file_name.endswith('_AF.png'):
            if file_name not in ['airfoil_plot_AF.png', 'visual_file_rotated_AF.txt', 'airfoil_coordinates_offset_AF.txt']:
                os.remove(file_name) 
