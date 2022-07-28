import matplotlib.pyplot as plt
import numpy as np

# Calculates analytical solutions for position with the given increasing-rotation
# velocity field. Plots the final positions and saves values as .npy

# Set initial x/y coordinates of detectors
x = 0.5 + 0.25 * np.arange(0, 100.0) / 100.0
y = np.zeros(100) + 0.5

# Set timestep parameters:
tmax = 8

r = 0.25 * np.arange(0, 100.0) / 100.0
theta = tmax**2.0 / 16.0
X = r * np.cos(theta) + 0.5
Y = r * np.sin(theta) + 0.5

plt.plot(X, Y, ".")
# print (x-X).max(), (y-Y).max()
plt.show()

# X.tofile('Xvals.txt',sep=' ')
# Y.tofile('Yvals.txt',sep=' ')
np.save("Xvals.npy", X)
np.save("Yvals.npy", Y)
