"""
Some testruns of the Data reg class
"""
import sys
sys.path.append("../")

import numpy as np
from random import random
from imageio import imread
import matplotlib.pyplot as plt
from Class.P1_Class import Data_reg


# Make data.
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
x, y = np.meshgrid(x, y)

def FrankeFunction(x, y):
    term1 = 0.75 * np.exp(-(0.25 * (9 * x - 2) ** 2) - 0.25 * ((9 * y - 2) ** 2))
    term2 = 0.75 * np.exp(-((9 * x + 1) ** 2) / 49.0 - 0.1 * (9 * y + 1))
    term3 = 0.5 * np.exp(-(9 * x - 7) ** 2 / 4.0 - 0.25 * ((9 * y - 3) ** 2))
    term4 = -0.2 * np.exp(-(9 * x - 4) ** 2 - (9 * y - 7) ** 2)
    return term1 + term2 + term3 + term4

z = FrankeFunction(x, y) + np.random.normal(0, 1, x.shape)

terrain = imread("../Data/n59_e010_1arc_v3.tif")
# Reducing the area we look at to make the analysis manageable
z = terrain[0:int(terrain.shape[1] / 5), 0:int(terrain.shape[1] / 5)]

# Show data image
plt.figure()
plt.title("Terrain of Finnemarka")
plt.imshow(z, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")

x = np.linspace(0, z.shape[1], z.shape[1])
y = np.linspace(0, z.shape[0], z.shape[0])
x, y = np.meshgrid(x, y)


A = Data_reg(x, y, z) # Initialising the class
A.Set_up_data(deg = 20)
A.PlaneOLSReg()
z_pred = A.Pred_data(A.X)

plt.figure()
plt.title("Pred of Finnemarka")
plt.imshow(z_pred, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")

plt.show()