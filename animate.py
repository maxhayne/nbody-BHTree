import numpy as np 
from math import sqrt
from matplotlib import pyplot as plt 
from matplotlib import animation

file_name = "nbody.mp4"

with open("locations.data") as file:
    bodies = int(file.readline())
    dimension = int(file.readline())
    dt = float(file.readline())
    duration = float(file.readline())
    mins = file.readline().split(',')
    xMin = float(mins[0])
    yMin = float(mins[1])
    sideLength = float(file.readline().strip())

    masses = []
    for i in range(bodies):
        masses.append(float(file.readline()))

     # frames per second for the video to be saved in
    fps = 60

    # total number of frames in the video
    total_frames = int(duration/dt)

    fig = plt.figure()
    sub_plot = fig.add_subplot(111, aspect='equal', autoscale_on=False,
    xlim=(xMin, xMin+sideLength), ylim=(yMin, yMin+sideLength))
    sub_plot.set_facecolor('black');
    fig.tight_layout(pad=0);

    # The positions of the bodies will be marked by blue circles of size 6
    positions, = sub_plot.plot([], [], 'wo', ms=1, alpha=0.50)

    # Trails will add red points at all the previous locations of the bodies
    trails, = sub_plot.plot([], [], 'r.', ms=2, alpha=0.20)

    # Store the history of the bodies
    x_history = []
    y_history = []

    def animate(i):
        """Animate function required by matplotlib.
        i - does nothing, but required by matplotlib"""

        xPosition = []
        yPosition = []

        # Read and store the masses of the bodies
        m = 0
        for point in file:
            x = float(point.split('\t')[0])
            y = float(point.split('\t')[1])
            #x_history.append(x)
            #y_history.append(y)
            xPosition.append(x)
            yPosition.append(y)
            m += 1
            if (m >= bodies): break

        #trails.set_data(x_history, y_history)
        positions.set_data(xPosition, yPosition)

        return positions#, trails

    anim = animation.FuncAnimation(fig, animate, frames=total_frames)
    anim.save(file_name, fps=fps, extra_args=['-vcodec', 'libx264'], dpi=300)