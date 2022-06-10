import numpy as np
import matplotlib.pyplot as plt

w = 1
w0 = np.linspace(0,2, 101)
f0 = 1
b = 0.1 * w
denom = np.square(w0 * w0 - w * w) + 4 * b * b * w * w
A2 = f0*f0 / denom
plt.plot(w0, A2, label='FWHM={}, Q={}'.format(2*b, w/(2*b)))

plt.legend()
plt.show()