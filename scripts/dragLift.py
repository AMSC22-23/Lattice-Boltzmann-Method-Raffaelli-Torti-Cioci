import matplotlib.pyplot as plt

# path to data file
file_path = 'outputs/lift_drag_out.txt'

frames = []
drags = []
lifts = []
# read while end of file is not reached
with open(file_path, 'r') as f:
    while True:
        line = f.readline()
        if not line:
            break
        frames.append(int(line.split()[0]))
        line = f.readline()
        drag = float(line.split()[0])
        lift = float(line.split()[1])
        drags.append(drag)
        lifts.append(lift)
        
# plot drag and lift in the same plot
plt.plot(frames, drags, label='drag')
plt.plot(frames, lifts, label='lift')
plt.xlabel('frame')
plt.legend()
plt.savefig('outputs/lift_drag.png')