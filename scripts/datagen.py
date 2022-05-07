from math import sin, cos
import random
import numpy as np


n = 2**7
X = np.linspace(0, 10, n)
Y = np.sin(X)

data_file = open("../data/generated_data.csv", "w")

for y in Y:
    data_file.write(f"{y:.4f}\n")

data_file.close()
