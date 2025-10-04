import numpy as np
from matplotlib import pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_no_Noise/Trajectory_i_0_5')



#read files: 
time = np.genfromtxt("time.txt")

MC_mean = np.genfromtxt("mean.txt")

i = 0.5
T = 10


############################     LR CIRCUIT     #########################################

def analytical_Voltage_mean(i):
    return np.where(i <= 1, 0, np.sqrt(i**2 -1))


fig, ax = plt.subplots(figsize = (6, 4.5))

ax.plot(time, MC_mean/(np.arcsin(i)), color = 'red', label = 'numerical')
#ax.plot(Relative_Current, analytical_Voltage_mean(Relative_Current), color = 'blue', label = 'analytical')
ax.set_title(r'Trajectory for $i = 0.5$')
ax.set_xlabel(r'$t /(t_0)$')
ax.set_ylabel(r'$\varphi /(\arcsin(i) \phi_0) $ ')
ax.grid(True)
ax.legend()



plt.tight_layout()
plt.savefig('JJ_R_Circuit_no_Noise_Trajectory_i_0_5.pdf')
plt.show()


