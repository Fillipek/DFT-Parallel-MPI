from matplotlib import pyplot as plt
import numpy as np

DATA_DIR = "data/"
FILENAME = "fft_data.csv"

file = open(DATA_DIR + FILENAME, "r")
X, Y_re, Y_im = [], [], []
for line in file.readlines():
    data_point = [float(x) for x in line.split(sep=",")]
    X.append(data_point[0])
    Y_re.append(data_point[1])
    Y_im.append(data_point[2])

Y_re = np.abs(Y_re)
Y_im = np.abs(Y_im)

y_max = np.max(np.concatenate((Y_re, Y_im), axis=0))
noise_thresh = 0.5 * y_max

plt.plot(X,Y_re)
plt.plot(X,Y_im)
plt.axhline(noise_thresh, color="r", linewidth=1, linestyle="dashed")
plt.show()

file.close()
