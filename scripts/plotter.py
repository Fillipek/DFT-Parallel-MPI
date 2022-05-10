from matplotlib import pyplot as plt
import numpy as np
import sys


files_to_plot = sys.argv[1:]

fig, axs = plt.subplots(1, len(files_to_plot))
plot_idx = 0;

for filename in files_to_plot:
    file = open(filename, "r")
    X, Y_re, Y_im = [], [], []
    for line in file.readlines():
        data_point = [float(x) for x in line.split(sep=",")]
        X.append(data_point[0])
        Y_re.append(data_point[1])
        Y_im.append(data_point[2])

    axs[plot_idx].set_title(filename)
    axs[plot_idx].plot(X,Y_re)
    axs[plot_idx].plot(X,Y_im)
    

    file.close()
    plot_idx += 1

plt.show()