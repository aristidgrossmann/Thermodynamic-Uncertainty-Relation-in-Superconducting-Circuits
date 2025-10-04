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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/Vortrag/figs/04/Ratchet')

sys.path.append( os.getcwd())

from Exact_Solutions import voltage_ratchet,  Uncertainty_Product_Ratchet, Conjecture_Ratchet


import Exact_Solutions 

Exact_Solutions.L = 1
Exact_Solutions.V_0 = 1
Exact_Solutions.a = Exact_Solutions.L*0.5
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a



i_0_min = -2
i_0_max = 2
N_i_0 = 100
i_0_vals = np.logspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.1, 0.3, 0.5])
D_labels = [r"$D = 0.1$", r"$D = 0.3$", r"$D = 0.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]




fig, ax = plt.subplots(figsize=(4.5, 3))

for i in range(len(D_vals)):
    Uncertainty_product_vals = np.ones_like(i_0_vals)
    Conjecture_vals = np.ones_like(i_0_vals)

    for j in range(len(Uncertainty_product_vals)): 
        Uncertainty_product_vals[j] = Uncertainty_Product_Ratchet(i_0 = i_0_vals[j], D = D_vals[i])
        Conjecture_vals[j] = Conjecture_Ratchet(i_0 = i_0_vals[j], D = D_vals[i])


    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    ax.plot(i_0_vals, Uncertainty_product_vals, linestyle = '-', label = D_labels[i], color = colors[i], zorder = 1)

    linestyles = ['dashdot', 'dashed', 'dotted']
    ax.plot(i_0_vals, Conjecture_vals, linestyle = 'dotted', color = colors[i], zorder = 1)




ax.hlines(2, 10**i_0_min, 10**i_0_max, linestyle = '--', linewidth = 2, color = 'black')



ax.set_xlabel('$i_0 / U_0$')

# ax.set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$', rotation=0)
# ax.yaxis.set_label_coords(-0.08, 0.5)
ax.set_yscale('log')
ax.set_xscale('log')

plt.tight_layout()

ax.set_xticks(
    ticks=[10**-2, 10**-1, 10**0, 10**1, 10**2],
    labels=[r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$']
)

ax.set_yticks(
    ticks=[2, 10],
    labels=[r'$2$', r'$10$']
)

from matplotlib.ticker import LogLocator
ax.xaxis.set_minor_locator(LogLocator(subs='all'))  # Auto minor ticks for log scale
ax.yaxis.set_minor_locator(LogLocator(subs='all'))
ax.tick_params(axis='both', which='minor', labelleft=False, labelbottom=False)  # Hide labels


ax.text(0.58, 0.97, r'$\tilde{\mathcal{Q}}$', transform=ax.transAxes,
         fontsize=11, va='top', ha='left')

ax.text(0.53, 0.65, r'$\tilde{\mathcal{L}}$', transform=ax.transAxes,
         fontsize=11, va='top', ha='left')

ax.legend(loc = 'upper left')

plt.savefig('TUR_Ratchet_Symmetric.svg')
plt.show()



