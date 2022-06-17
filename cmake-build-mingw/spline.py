import matplotlib.pyplot as plt
import numpy as np

dop_array = open("spline1.txt").read().split()
spline1_x = np.array([float(dop_array[x]) for x in range(0, len(dop_array), 2)])
spline1_y = np.array([float(dop_array[y]) for y in range(1, len(dop_array), 2)])
dop_array = open("spline2.txt").read().split()
spline2_x = np.array([float(dop_array[x]) for x in range(0, len(dop_array), 2)])
spline2_y = np.array([float(dop_array[y]) for y in range(1, len(dop_array), 2)])

plt.plot(spline1_x, spline1_y, color="red", label="Spline1")
plt.plot(spline2_x, spline2_y, color="blue", label="Spline2")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()



