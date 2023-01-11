import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
from matplotlib import cm


matplotlib.use('TkAgg')
plt.interactive(False)

def parse_array(file_path: str):

    X = []
    Y = []
    Z = []
    density = []

    with open(file_path) as map_file:

        for line in map_file:
            values = line.split(" ")

            x = float(values[0])
            y = float(values[1])
            z = float(values[2])
            dens = float(values[3])

            if dens < 0:
                continue

            X.append(x)
            Y.append(y)
            Z.append(z)
            density.append(dens)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    img = ax.scatter(X, Y, Z, c=density, s=density)
    plt.colorbar(img)
    plt.show()


parse_array("./debug/aligned_fragments/array3d.csv")
