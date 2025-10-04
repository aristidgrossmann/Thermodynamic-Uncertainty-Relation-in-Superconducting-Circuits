import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os
import sys

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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/JJ')

sys.path.append( os.getcwd())

from Exact_Solutions import L, Mean_Velocity, Variance_Velocity, Uncertainty_product, Conjecture




i_0_numerical = np.genfromtxt('i_0_vals_numerical.txt')

i_0_min = 0.01
i_0_max = np.max(i_0_numerical)
N_i_0 = 100
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.1, 0.5, 1.5])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]




fig, ax = plt.subplots(figsize=(4.5, 3))

for i in range(len(D_vals)):
    Uncertainty_product_vals = np.ones_like(i_0_vals)
    Conjecture_vals = np.ones_like(i_0_vals)

    for j in range(len(Uncertainty_product_vals)): 
        Uncertainty_product_vals[j] = Uncertainty_product(i_0 = i_0_vals[j], D = D_vals[i])
        Conjecture_vals[j] = Conjecture(i_0 = i_0_vals[j], D = D_vals[i])
        print(j)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    ax.plot(i_0_vals, Uncertainty_product_vals, linestyle = '-', color = colors[i], zorder = 1)

    linestyles = ['dashdot', 'dashed', 'dotted']
    ax.plot(i_0_vals, Conjecture_vals, linestyle = 'dotted', color = colors[i], zorder = 1)


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/JJ')

Numerical_mean_velocitites_D_0_1 = np.genfromtxt('0.1_mean_velocities.txt')
# for i in range(len(i_0_numerical)):
#     Numerical_mean_velocitites_D_0_1[i] = Mean_Velocity(i_0 = i_0_numerical[i], D = D_vals[0])

Numerical_mean_velocitites_D_0_5 = np.genfromtxt('0.5_mean_velocities.txt')
Numerical_mean_velocitites_D_1_5 = np.genfromtxt('1.5_mean_velocities.txt')

Numerical_variance_velocitites_D_0_1 = np.genfromtxt('0.1_variance_velocities.txt')

# for i in range(len(i_0_numerical)):
#     Numerical_variance_velocitites_D_0_1[i] = Variance_Velocity(i_0 = i_0_numerical[i], D = D_vals[0])*1/(1+1/(i_0_vals[i]**4 + 10**3))

Numerical_variance_velocitites_D_0_5 = np.genfromtxt('0.5_variance_velocities.txt')
Numerical_variance_velocitites_D_1_5 = np.genfromtxt('1.5_variance_velocities.txt')

TUR_vals_numerical_D_0_1 = Numerical_variance_velocitites_D_0_1/Numerical_mean_velocitites_D_0_1*i_0_numerical/D_vals[0]
TUR_vals_numerical_D_0_5 = Numerical_variance_velocitites_D_0_5/Numerical_mean_velocitites_D_0_5*i_0_numerical/D_vals[1]
TUR_vals_numerical_D_1_5 = Numerical_variance_velocitites_D_1_5/Numerical_mean_velocitites_D_1_5*i_0_numerical/D_vals[2]



ax.scatter(i_0_numerical[::2], TUR_vals_numerical_D_0_1[::2], marker = 'o' , s = 20, color = colors[0], zorder = 2)
ax.scatter(i_0_numerical[::2], TUR_vals_numerical_D_0_5[::2], marker = 'o' , s = 20, color = colors[1], zorder = 2)
ax.scatter(i_0_numerical[::2], TUR_vals_numerical_D_1_5[::2], marker = 'o' , s = 20, color = colors[2], zorder = 2)



ax.hlines(2, i_0_min - 0.05, i_0_max+0.05, linestyle = '--', linewidth = 2, color = 'black')



ax.set_xlabel('$I_0 / I_{\mathrm{c}}$')
ax.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
ax.yaxis.set_label_coords(-0.08, 0.7)
ax.set_yscale('log')

plt.tight_layout()
plt.legend()

ax.set_xticks(
    ticks=[0, 1, 2, 3],
    labels=[r'$0$', r'$1$', r'$2$', r'$3$']
)

ax.set_yticks(
    ticks=[2, 10],
    labels=[r'$2$', r'$10^1$']
)



plt.savefig('JJ_Uncertainty_Product.pdf')
plt.show()



