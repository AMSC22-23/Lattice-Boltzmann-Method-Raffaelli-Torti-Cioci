import math

def generate_circle_coordinates(center_x, center_y, radius, grid_size):
    circle_coordinates = []
    for x in range(grid_size[0]):
        for y in range(grid_size[1]):
            # Calculate the distance from the center to the point
            distance = math.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
            
            # If the distance is less than or equal to the radius, the point is inside the circle
            if distance <= radius:
                circle_coordinates.append((x, y))

    return circle_coordinates

# Example usage for a 200x50 grid
grid_size = (200, 50)
center = (100, 25)
radius = 12

circle_coordinates = generate_circle_coordinates(center[0], center[1], radius, grid_size)

# print coordinates to file separated by spaces and \n
with open('circle_coordinates.txt', 'w') as f:
    for x, y in circle_coordinates:
        f.write(f'{x} {y}\n')