import numpy as np
from matplotlib import pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/LR_Circuit')



#read files: 
time = np.genfromtxt("time.txt")

MC_mean = np.genfromtxt("mean.txt")
MC_mean_std = np.genfromtxt("mean_std.txt")

MC_variance = np.genfromtxt("variance.txt")
MC_variance_std = np.genfromtxt("variance_std.txt")


############################     LR CIRCUIT     #########################################

def analytical_mean(t):
    return np.zeros_like(t)

def analytical_variance(t):
    return 0.5*(1- np.exp(-2*t))



fig, (ax1, ax2) = plt.subplots(2,1, figsize = (6, 6))

ax1.plot(time, MC_mean, color = 'red', label = 'numerical')
ax1.fill_between(time, MC_mean - MC_mean_std, MC_mean + MC_mean_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax1.plot(time, analytical_mean(time), color = 'blue', label = 'analytical')
ax1.set_title('Mean')
ax1.set_xlabel(r'$t \cdot R/L$')
ax1.set_ylabel(r'$\langle \phi \rangle / \sqrt{2 k_{\mathrm{B}}T L}$ ')
ax1.grid(True)
ax1.legend()

ax2.plot(time, MC_variance, color = 'red', label = 'numerical')
ax2.fill_between(time, MC_variance - MC_variance_std, MC_variance + MC_variance_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax2.plot(time, analytical_variance(time), color = 'blue', label = 'analytical')
ax2.set_title('Variance')
ax2.set_xlabel(r'$t \cdot R/L$')
ax2.set_ylabel(r'$\langle \langle \phi^2 \rangle \rangle / (2 k_{\mathrm{B}}T L)$')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig('result_LC_Circuit.pdf')
plt.show()


