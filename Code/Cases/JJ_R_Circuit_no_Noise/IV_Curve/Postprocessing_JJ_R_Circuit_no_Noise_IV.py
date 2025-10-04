import numpy as np
from matplotlib import pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_no_Noise/IV_Curve')



#read files: 
Relative_Current = np.genfromtxt("Relative_Current_IV.txt")

Voltage_mean = np.genfromtxt("Voltage_mean_IV.txt")
Voltage_mean_std = np.genfromtxt("Voltage_mean_std_IV.txt")



############################     LR CIRCUIT     #########################################

def analytical_Voltage_mean(i):
    return np.where(i <= 1, 0, np.sqrt(i**2 -1))


fig, ax = plt.subplots(figsize = (6, 4.5))

ax.plot(Relative_Current, Voltage_mean, color = 'red', label = 'numerical')
ax.fill_between(Relative_Current, Voltage_mean - Voltage_mean_std, Voltage_mean + Voltage_mean_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax.plot(Relative_Current, analytical_Voltage_mean(Relative_Current), color = 'blue', label = 'analytical')
ax.set_title('I-V curve')
ax.set_xlabel(r'$I_0 / I_{\mathrm{c}}$')
ax.set_ylabel(r'$V/(R I_{\mathrm{c}})$ ')
ax.grid(True)
ax.legend()



plt.tight_layout()
plt.savefig('JJ_R_Circuit_no_Noise_IV.pdf')
plt.show()


