import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/TUR')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, T_1_paper, Delta_T_2_paper, Delta_T_2_me


i_0_min = 0.1
i_0_max = 2
N_i_0 = 58
N_i_0_numerical = 20
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)
i_0_numerical = np.linspace(i_0_min, i_0_max, N_i_0_numerical)

D_vals = np.array([0.1, 0.5, 1, 2, 5, 10])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1$", r"$D = 2$", r"$D = 5$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))




fig, ax = plt.subplots(figsize = (7, 5.2))

for i in range(len(D_vals)):
    TUR_vals = np.ones_like(i_0_vals)
    for j in range(len(TUR_vals)): 
        TUR_vals[j] = L*Delta_T_2_paper(i_0 = i_0_vals[j], D = D_vals[i]) / T_1_paper(i_0 = i_0_vals[j], D = D_vals[i])**2 * i_0_vals[j]/D_vals[i]


    ax.plot(i_0_vals, TUR_vals, label=D_labels[i] + ' analytical', color = colors[i+1])


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/Components')

Numerical_mean_velocitites_D_0_1 = np.genfromtxt('0.1_mean_velocities.txt')
Numerical_mean_velocitites_D_0_5 = np.genfromtxt('0.5_mean_velocities.txt')
Numerical_mean_velocitites_D_1 = np.genfromtxt('1.0_mean_velocities.txt')
Numerical_mean_velocitites_D_2 = np.genfromtxt('2.0_mean_velocities.txt')
Numerical_mean_velocitites_D_5 = np.genfromtxt('5.0_mean_velocities.txt')
Numerical_mean_velocitites_D_10 = np.genfromtxt('10._mean_velocities.txt')

Numerical_variance_velocitites_D_0_1 = np.genfromtxt('0.1_variance_velocities.txt')
Numerical_variance_velocitites_D_0_5 = np.genfromtxt('0.5_variance_velocities.txt')
Numerical_variance_velocitites_D_1 = np.genfromtxt('1.0_variance_velocities.txt')
Numerical_variance_velocitites_D_2 = np.genfromtxt('2.0_variance_velocities.txt')
Numerical_variance_velocitites_D_5 = np.genfromtxt('5.0_variance_velocities.txt')
Numerical_variance_velocitites_D_10 = np.genfromtxt('10._variance_velocities.txt')

TUR_vals_numerical_D_0_1 = Numerical_variance_velocitites_D_0_1/Numerical_mean_velocitites_D_0_1*i_0_numerical/D_vals[0]
TUR_vals_numerical_D_0_5 = Numerical_variance_velocitites_D_0_5/Numerical_mean_velocitites_D_0_5*i_0_numerical/D_vals[1]
TUR_vals_numerical_D_1 = Numerical_variance_velocitites_D_1/Numerical_mean_velocitites_D_1*i_0_numerical/D_vals[2]
TUR_vals_numerical_D_2 = Numerical_variance_velocitites_D_2/Numerical_mean_velocitites_D_2*i_0_numerical/D_vals[3]
TUR_vals_numerical_D_5 = Numerical_variance_velocitites_D_5/Numerical_mean_velocitites_D_5*i_0_numerical/D_vals[4]
TUR_vals_numerical_D_10 = Numerical_variance_velocitites_D_10/Numerical_mean_velocitites_D_10*i_0_numerical/D_vals[5]



ax.scatter(i_0_numerical, TUR_vals_numerical_D_0_1, marker = 'x' , s = 25, label=D_labels[0] + ' numerical', color = 'red')
ax.scatter(i_0_numerical, TUR_vals_numerical_D_0_5, marker = '+' , s = 50, label=D_labels[1] + ' numerical', color = 'red')
ax.scatter(i_0_numerical, TUR_vals_numerical_D_1,  marker = '^' , s = 25, label=D_labels[2] + ' numerical', color = 'red')
ax.scatter(i_0_numerical, TUR_vals_numerical_D_2, marker = 'o' , s = 25, label=D_labels[3] + ' numerical', color = 'red')
ax.scatter(i_0_numerical, TUR_vals_numerical_D_5, marker = '<' , s = 25, label=D_labels[4] + ' numerical', color = 'red')
ax.scatter(i_0_numerical, TUR_vals_numerical_D_10, marker = '>' , s = 25, label=D_labels[5] + ' numerical', color = 'red')



ax.hlines(2, i_0_min - 0.05, i_0_max+0.05, linestyle = '--', linewidth = 2, label = 'TUR lower bound', color = 'black')



ax.set_xlabel(r"$i_0$", fontsize=12)
ax.set_ylabel("Uncertainty Product", fontsize=12)
ax.grid(True)
ax.set_yscale('log')

plt.tight_layout()
plt.legend()


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/TUR')

plt.savefig('TUR.pdf')
plt.show()



