import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import matplotlib as mpl


mpl.rcParams.update({
    "text.usetex": False,            # Use matplotlib's internal math rendering
    "mathtext.fontset": "cm",        # Use Computer Modern (LaTeX-like)
    "font.family": "serif",          # Serif font family
    "font.size": 11,                 # Base font size
    "xtick.labelsize": 10,           # Tick label font size
    "ytick.labelsize": 10,
    "axes.labelsize": 11,            # Axis label font size
    "legend.fontsize": 10,           # Legend font size (default: 10)
})

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, T_1_paper, Delta_T_2_paper, Delta_T_2_me




os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/Components')

i_0_numerical = np.genfromtxt('i_0_vals_numerical.txt')



i_0_vals = np.linspace(np.min(i_0_numerical), np.max(i_0_numerical), 200)

D_vals = np.array([0.1, 0.5, 1.5])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = ["black", "blue", "red"]




fig, ax = plt.subplots(figsize = (4.5, 3))

for i in range(len(D_vals)):
    mean_vals = np.ones_like(i_0_vals)
    for j in range(len(mean_vals)):
        mean_vals[j] = L /T_1_paper(i_0 = i_0_vals[j], D = D_vals[i])

    ax.plot(i_0_vals, mean_vals, label=D_labels[i] + ' analytical', color = colors[i])



Numerical_mean_velocitites_D_0_1 = np.genfromtxt('0.1_mean_velocities.txt')
Numerical_mean_velocitites_D_0_5 = np.genfromtxt('0.5_mean_velocities.txt')
Numerical_mean_velocitites_D_1 = np.genfromtxt('1.5_mean_velocities.txt')


increments = 2
ax.scatter(i_0_numerical[::increments], Numerical_mean_velocitites_D_0_1[::increments], marker = 'o' , s = 10, label=D_labels[0] + ' numerical', color = colors[0])
ax.scatter(i_0_numerical[::increments], Numerical_mean_velocitites_D_0_5[::increments], marker = 'o' , s = 10, label=D_labels[1] + ' numerical', color = colors[1])
ax.scatter(i_0_numerical[::increments], Numerical_mean_velocitites_D_1[::increments],  marker = 'o' , s = 10, label=D_labels[2] + ' numerical', color = colors[2])


def w_0(i_0):
    if i_0 < 1:
        return 0
    else: 
        return np.sqrt(i_0**2 -1)
    
no_noise = np.ones_like(i_0_vals)
for i in range(len(i_0_vals)):
    no_noise[i] = w_0(i_0_vals[i])

ax.plot(i_0_vals, no_noise, linestyle = 'dotted',  color = "purple")
ax.plot(i_0_vals, i_0_vals, linestyle = 'dashed',  color = "purple")
    

ax.set_xlabel(r"$i_0$")
ax.set_ylabel(r'$\langle v\rangle$')

ax.set_xticks([0,1,2,3])
ax.set_yticks([0,1,2,3])

plt.tight_layout()

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components')

plt.savefig('Mean_Velocity.pdf')
plt.show()




