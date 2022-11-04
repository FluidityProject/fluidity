import matplotlib.pyplot as plt
import numpy as np

num_detectors = 100
x = 0.5 + 0.25 * np.arange(0, float(num_detectors)) / float(num_detectors)
y = np.zeros(num_detectors) + 0.5

t = 0.0
n_cycles = 1
dt = 0.1 / n_cycles
tmax = 8


def vel(x, y):
    return [-(y - 0.5), x - 0.5]


while t < tmax:
    t = t + dt
    [k1_x, k1_y] = vel(x, y)
    [k2_x, k2_y] = vel(x + 0.5 * dt * k1_x, y + 0.5 * dt * k1_y)
    [k3_x, k3_y] = vel(x + 0.5 * dt * k2_x, y + 0.5 * dt * k2_y)
    [k4_x, k4_y] = vel(x + dt * k3_x, y + dt * k3_y)
    x = x + dt * (k1_x / 6.0 + k2_x / 3.0 + k3_x / 3.0 + k4_x / 6.0)
    y = y + dt * (k1_y / 6.0 + k2_y / 3.0 + k3_y / 3.0 + k4_y / 6.0)

plt.plot(x, y, ".")
# show()

x.tofile("Xvals.txt", sep=" ")
y.tofile("Yvals.txt", sep=" ")
