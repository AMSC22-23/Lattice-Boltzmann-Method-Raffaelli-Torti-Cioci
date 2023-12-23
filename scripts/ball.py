import math

def generate_circle_coordinates(center_x, center_y, radius, grid_size):
    circle_coordinates = []
    for theta in range(0, 360):
        # Convert degrees to radians
        radian_theta = math.radians(theta)
        
        # Calculate coordinates
        x = center_x + int(radius * math.cos(radian_theta))
        y = center_y + int(radius * math.sin(radian_theta))
        
        # Ensure the coordinates are within the grid boundaries
        if 0 <= x < grid_size[0] and 0 <= y < grid_size[1]:
            circle_coordinates.append((x, y))

    return circle_coordinates

# Example usage for a 200x50 grid
grid_size = (200, 50)
center = (20, 25)
radius = 12

circle_coordinates = generate_circle_coordinates(center[0], center[1], radius, grid_size)

# print coordinates to file separated by spaces and \n
with open('circle.txt', 'w') as f:
    for x, y in circle_coordinates:
        f.write(f'{x} {y}\n')