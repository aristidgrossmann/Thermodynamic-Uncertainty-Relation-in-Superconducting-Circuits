import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components_SDE_Solver')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, T_1_paper, Delta_T_2_paper, Delta_T_2_me


i_0_min = 0.1
i_0_max = 2
N_i_0 = 58
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.1, 0.5, 1, 2, 5, 10])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1$", r"$D = 2$", r"$D = 5$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))




fig, ax = plt.subplots(figsize = (7, 5.2))

for i in range(len(D_vals)):
    Variance_vals = np.ones_like(i_0_vals)
    for j in range(len(Variance_vals)):
        Variance_vals[j] = L**2*Delta_T_2_paper(i_0 = i_0_vals[j], D = D_vals[i]) / T_1_paper(i_0 = i_0_vals[j], D = D_vals[i])**3

    ax.plot(i_0_vals, Variance_vals, label=D_labels[i] + ' analytical', color = colors[i+1])


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/Components_SDE')

i_0_numerical = np.genfromtxt('i_0_vals_numerical.txt')

Numerical_variance_velocitites_D_0_1 = np.genfromtxt('0.1_variance_velocities.txt')
Numerical_variance_velocitites_std_D_0_1 = np.genfromtxt('0.1_variance_velocities_std.txt')

Numerical_variance_velocitites_D_0_5 = np.genfromtxt('0.5_variance_velocities.txt')
Numerical_variance_velocitites_std_D_0_5 = np.genfromtxt('0.5_variance_velocities_std.txt')

Numerical_variance_velocitites_D_1 = np.genfromtxt('1.0_variance_velocities.txt')
Numerical_variance_velocitites_std_D_1 = np.genfromtxt('1.0_variance_velocities_std.txt')

Numerical_variance_velocitites_D_2 = np.genfromtxt('2.0_variance_velocities.txt')
Numerical_variance_velocitites_std_D_2 = np.genfromtxt('2.0_variance_velocities_std.txt')

Numerical_variance_velocitites_D_5 = np.genfromtxt('5.0_variance_velocities.txt')
Numerical_variance_velocitites_std_D_5 = np.genfromtxt('5.0_variance_velocities_std.txt')

Numerical_variance_velocitites_D_10 = np.genfromtxt('10._variance_velocities.txt')
Numerical_variance_velocitites_std_D_10 = np.genfromtxt('10._variance_velocities_std.txt')


ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_0_1, yerr = Numerical_variance_velocitites_std_D_0_1, 
            fmt = '.', capsize=5, label=D_labels[0] + ' numerical', color = 'red')

ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_0_5, yerr = Numerical_variance_velocitites_std_D_0_5, 
            fmt = '.', capsize=5, label=D_labels[1] + ' numerical', color = 'red')

ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_1, yerr = Numerical_variance_velocitites_std_D_1, 
            fmt = '.', capsize=5, label=D_labels[2] + ' numerical', color = 'red')

ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_2, yerr = Numerical_variance_velocitites_std_D_2, 
            fmt = '.', capsize=5, label=D_labels[3] + ' numerical', color = 'red')

ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_5, yerr = Numerical_variance_velocitites_std_D_5, 
            fmt = '.', capsize=5, label=D_labels[4] + ' numerical', color = 'red')

ax.errorbar(i_0_numerical, Numerical_variance_velocitites_D_10, yerr = Numerical_variance_velocitites_std_D_10, 
            fmt = '.', capsize=5, label=D_labels[5] + ' numerical', color = 'red')


ax.set_xlabel(r"$i_0$", fontsize=12)
ax.set_ylabel("Variance Velo", fontsize=12)
ax.grid(True)
plt.tight_layout()
plt.legend()

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components_SDE_Solver')

plt.savefig('Variance_Velocity_MC.pdf')
plt.show()



