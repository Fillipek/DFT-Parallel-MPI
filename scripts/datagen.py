from math import sin, cos
from random import random
import numpy as np
from matplotlib import pyplot as plt


n = 2**10
X = np.linspace(0, 12, n)
Y = np.sin(np.copy(X)) + np.cos(2 * np.copy(X)) + 0.3*np.random.random(n)

data_file = open("data/generated_data.csv", "w")

for y in Y:
    data_file.write(f"{y:.4f}\n")

data_file.close()

plt.title("Generated data");
plt.plot(X,Y)
plt.show()
