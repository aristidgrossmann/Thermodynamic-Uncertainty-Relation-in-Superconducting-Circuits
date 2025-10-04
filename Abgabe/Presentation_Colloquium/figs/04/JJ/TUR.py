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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/Vortrag/figs/04/JJ')

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
    ax.plot(i_0_vals, Uncertainty_product_vals, linestyle = '-', label = D_labels[i], color = colors[i], zorder = 1)



ax.hlines(2, i_0_min - 0.05, i_0_max+0.05, linestyle = '--', linewidth = 2, color = 'black')


ax.text(0.32, 0.97, r'$\tilde{\mathcal{Q}}$', transform=ax.transAxes,
         fontsize=11, va='top', ha='left')

ax.set_xlabel('$I_0 / I_{\mathrm{c}}$')
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




plt.savefig('JJ_TUR.svg')
plt.show()



