# this script is a demo of FFT in Python.
# Can be used as a reference to check if C implementation is correct.

import numpy as np
from matplotlib import pyplot as plt

n = 2**7
X = np.linspace(0, 10, n)
Y = np.sin(np.copy(X))

FFT = np.fft.fft(Y, n)
IFFT = np.fft.ifft(FFT, n)

plt.subplot(1,3,1); plt.title("Input")
plt.plot(X, Y)

plt.subplot(1,3,2); plt.title("FFT")
plt.plot(X, FFT.real)
plt.plot(X, FFT.imag)

plt.subplot(1,3,3); plt.title("Inversed FFT")
plt.plot(X, IFFT.real)
plt.plot(X, IFFT.imag)

plt.show()
