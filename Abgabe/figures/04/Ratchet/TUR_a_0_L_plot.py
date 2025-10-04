import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import LogLocator

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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/Ratchet')

sys.path.append( os.getcwd())

from Exact_Solutions import  Uncertainty_Product_Ratchet_a_0, Conjecture_Ratchet_a_0, Uncertainty_Product_Ratchet_a_L, Conjecture_Ratchet_a_L


import Exact_Solutions 

Exact_Solutions.L = 1
Exact_Solutions.V_0 = 1
Exact_Solutions.a = Exact_Solutions.L*1.0
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a



i_0_min = -2
i_0_max = 1.82
N_i_0 = 100
i_0_vals = np.logspace(i_0_min, i_0_max, N_i_0)

i_0_vals = np.genfromtxt("i_0_vals.txt")

D_vals = np.array([0.1, 0.3, 0.5])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]




# fig, (ax1, ax2) = plt.subplots(2, 1,   figsize=(4.5, 5.4), sharex=True)
fig, (ax1, ax2) = plt.subplots(2, 1,   figsize=(4.5, 4.5), sharex=True)

for i in range(len(D_vals)):
    filename = "Uncertainty_Product_a_0_" + str(D_vals[i]) + ".txt"
    Uncertainty_product_vals = np.genfromtxt(filename)

    filename = "Conjecture_a_0_" + str(D_vals[i]) + ".txt"
    Conjecture_vals = np.genfromtxt(filename)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    ax1.plot(i_0_vals, Uncertainty_product_vals, linestyle = '-', color = colors[i], zorder = 1)

    linestyles = ['dashdot', 'dashed', 'dotted']
    ax1.plot(i_0_vals, Conjecture_vals, linestyle = 'dotted', color = colors[i], zorder = 1)



for i in range(len(D_vals)):

    filename = "Uncertainty_Product_a_L_" + str(D_vals[i]) + ".txt"
    Uncertainty_product_vals = np.genfromtxt(filename)

    filename = "Conjecture_a_L_" + str(D_vals[i]) + ".txt"
    Conjecture_vals = np.genfromtxt(filename)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    ax2.plot(i_0_vals, Uncertainty_product_vals, linestyle = '-', color = colors[i], zorder = 1)

    linestyles = ['dashdot', 'dashed', 'dotted']
    ax2.plot(i_0_vals, Conjecture_vals, linestyle = 'dotted', color = colors[i], zorder = 1)





ax1.hlines(2, 10**i_0_min, 10**i_0_max, linestyle = '--', linewidth = 2, color = 'black')
ax2.hlines(2, 10**i_0_min, 10**i_0_max, linestyle = '--', linewidth = 2, color = 'black')


# ax1.set_xlabel('$i_0 / U_0$')
ax1.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
# ax1.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$', rotation=0)
# ax1.yaxis.set_label_coords(-0.15, 0.5)
ax1.set_yscale('log')
ax1.set_xscale('log')


ax2.set_xlabel('$i_0 / U_0$')
ax2.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
# ax2.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$', rotation=0)
# ax2.yaxis.set_label_coords(-0.15, 0.5)
ax2.set_yscale('log')
ax2.set_xscale('log')


ax1.set_xticks(
    ticks=[10**-2, 10**-1, 10**0, 10**1, 10**2],
    labels=[r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$']
)

ax1.set_yticks(
    ticks=[2, 10**1, 10**2, 10**3],
    labels=[r'$2$', r'$10^1$', r'$10^2$', r'$10^3$']
)


ax2.set_xticks(
    ticks=[10**-2, 10**-1, 10**0, 10**1, 10**2],
    labels=[r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$']
)

ax2.set_yticks(
    ticks=[2, 10],
    labels=[r'$2$', r'$10$']
)

ax1.xaxis.set_minor_locator(LogLocator(subs='all'))  # Auto minor ticks for log scale
ax1.yaxis.set_minor_locator(LogLocator(subs='all'))
ax1.tick_params(axis='both', which='minor', labelleft=False, labelbottom=False)  # Hide labels

ax2.xaxis.set_minor_locator(LogLocator(subs='all'))  # Auto minor ticks for log scale
ax2.yaxis.set_minor_locator(LogLocator(subs='all'))
ax2.tick_params(axis='both', which='minor', labelleft=False, labelbottom=False)  # Hide labels


ax1.text(0.03, 0.97, r'(a): $a = 0$', transform=ax1.transAxes,
         fontsize=11, va='top', ha='left')

ax2.text(0.03, 0.97, r'(b): $a = L$', transform=ax2.transAxes,
         fontsize=11, va='top', ha='left')

plt.tight_layout()


plt.savefig('TUR_Ratchet_a_0_L.pdf')
plt.show()



