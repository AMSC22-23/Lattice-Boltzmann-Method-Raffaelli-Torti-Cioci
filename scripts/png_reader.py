from PIL import Image

def find_white_pixels(image_path):
    # Open the image file
    img = Image.open(image_path)
    # Convert the image to black and white
    img = img.convert('1')
    # Get the size of the image
    width, height = img.size
    # Initialize an empty list to hold the coordinates of white pixels
    white_pixels = []

    # Loop over every pixel in the image
    for y in range(height):
        for x in range(width):
            # If the pixel is white (1 in a '1' mode image), add its coordinates to the list
            if img.getpixel((x, y)):
                white_pixels.append((x, y))

    return white_pixels

def write_coordinates_to_file(coordinates, output_file):
    # Open the output file in write mode
    with open(output_file, 'w') as f:
        # Write each coordinate to the file
        for coord in coordinates:
            f.write(f'{coord[0]} {coord[1]}\n')

# Use the functions
white_pixels = find_white_pixels('outputs/me.png')
write_coordinates_to_file(white_pixels, 'outputs/me.txt')