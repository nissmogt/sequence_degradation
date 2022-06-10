import numpy as np
import matplotlib.pyplot as plt

mu = 0
sigma = 0.4

x_all = np.arange(-1, 1, 0.01) # entire range of x, both in and out of spec
# mean = 0, stddev = 1, since Z-transform was calculated
aa = np.square(x_all-mu) / (2*np.square(sigma))
bb = sigma*np.sqrt(2*np.pi)
y = np.exp(-aa)/bb
# build the plot
fig, ax = plt.subplots(figsize=(9,6))
plt.style.use('fivethirtyeight')
ax.plot(x_all,y, c='black')
# plt.show()

# ax.fill_between(x,y,0, alpha=0.3, color='b')
ax.fill_between(x_all,y,0, alpha=0.1)
# ax.set_xlim([-4,4])
# ax.set_xlabel('# of Standard Deviations Outside the Mean')
# ax.set_yticklabels([])
# ax.set_title('Normal Gaussian Curve')

plt.savefig('normal_curve.svg', dpi=200, bbox_inches='tight')
plt.show()