import numpy as np
from matplotlib import pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/\Trajectory_i_0_5_ratio_1')



#read files: 
time = np.genfromtxt("time.txt")

MC_mean = np.genfromtxt("mean.txt")
MC_mean_std = np.genfromtxt("mean_std.txt")

MC_variance = np.genfromtxt("variance.txt")
MC_variance_std = np.genfromtxt("variance_std.txt")

i = 0.5
T = 10


fig, (ax1, ax2) = plt.subplots(2,1, figsize = (6, 6))

ax1.plot(time, MC_mean, color = 'red', label = 'numerical')
ax1.fill_between(time, MC_mean - MC_mean_std, MC_mean + MC_mean_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax1.set_title(r'Mean for $i = 0.5$; noise ratio = 1')
ax1.set_xlabel(r'$t /(t_0)$')
ax1.set_ylabel(r'$\varphi $ ')
ax1.grid(True)
ax1.legend()

ax2.plot(time, MC_variance, color = 'red', label = 'numerical')
ax2.fill_between(time, MC_variance - MC_variance_std, MC_variance + MC_variance_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax2.set_title('Variance for $i = 0.5$; noise ratio = 1')
ax2.set_xlabel(r'$t/t_0$')
ax2.set_ylabel(r'$\langle \langle \varphi^2 \rangle \rangle$')
ax2.grid(True)
ax2.legend()


plt.tight_layout()
plt.savefig('JJ_R_Circuit_Noise_Trajectory_i_0_5_ratio_1.pdf')
plt.show()




