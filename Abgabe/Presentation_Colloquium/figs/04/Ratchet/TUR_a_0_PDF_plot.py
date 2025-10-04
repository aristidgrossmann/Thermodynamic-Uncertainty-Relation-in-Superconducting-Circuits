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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/Vortrag/figs/04/Ratchet')

sys.path.append( os.getcwd())

from Exact_Solutions import T_1_Ratchet_a_0, T_1_Ratchet_a_L, I_plus_minus_Ratchet_a_0, I_minus_plus_Ratchet_a_0, I_plus_minus_Ratchet_a_L, I_minus_plus_Ratchet_a_L


import Exact_Solutions 

Exact_Solutions.L = 1
Exact_Solutions.V_0 = 1




i_0 = 6
N_phi = 1000
phi_vals = np.linspace(-0.5, 0.5, N_phi)



# D_vals = np.array([0.1, 0.3, 0.5])
D_vals = np.array([0.1])

D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]




# fig, (ax1, ax2) = plt.subplots(2, 1,   figsize=(4.5, 5.4), sharex=True)
fig, axs = plt.subplots(2, 2,   figsize=(5.4, 4.5), sharex='col')

### ax11
Exact_Solutions.a = Exact_Solutions.L*0.0
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
Exact_Solutions.V = Exact_Solutions.V_Ratchet
eps = 1e-3
phi_V_vals = np.linspace(-0.5 - eps, 0.5 + eps, 1000)
V_vals = np.ones_like(phi_V_vals)
for j in range(len(V_vals)):
    V_vals[j] = Exact_Solutions.V(phi_V_vals[j], i_0, 1)
axs[0,0].plot(phi_V_vals, V_vals, label = r'$U_{-}$', linestyle = '-',  color = "black")
axs[0,0].plot(phi_V_vals, -V_vals -V_vals[0], label = r'$U_{+}$', linestyle = 'dotted',  color = "red")

axs[0,0].legend(loc = 'upper right')
### ax12
Exact_Solutions.a = Exact_Solutions.L*1.0
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
Exact_Solutions.V = Exact_Solutions.V_Ratchet
for j in range(len(V_vals)):
    V_vals[j] = Exact_Solutions.V(phi_V_vals[j], i_0, 1)
axs[0,1].plot(phi_V_vals, V_vals, label = r'$U_{-}$',  color = "black")
axs[0,1].plot(phi_V_vals, -V_vals - V_vals[0], label = r'$U_{+}$', linestyle = 'dotted',  color = "red")

# axs[0,1].legend(loc = 'upper right')

for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)

    filename = "PDF_minus_a_0_" + str(D_vals[i]) + ".txt"

    PDF_minus = np.genfromtxt(filename)

    ## normalize
    PDF_minus = PDF_minus*len(PDF_minus)/np.sum(PDF_minus)
        

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    axs[1,0].plot(phi_vals, PDF_minus, label = r'$\tilde{J}_{-}$', linestyle = '-', color = colors[i], zorder = 1)
    axs[1,0].plot(phi_vals, PDF_minus[::-1], label = r'$\tilde{J}_{+}$', linestyle = 'dotted', color = "red", zorder = 1)

axs[1,0].legend(loc = 'upper right')



for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)


    filename = "PDF_minus_a_L_" + str(D_vals[i]) + ".txt"
    PDF_minus = np.genfromtxt(filename)
        
    ## normalize
    PDF_minus = PDF_minus*len(PDF_minus)/np.sum(PDF_minus)/(2*np.pi)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    axs[1,1].plot(phi_vals, PDF_minus, label = r'$\tilde{J}_{-}$',  linestyle = '-', color = colors[i], zorder = 1)
    axs[1,1].plot(phi_vals, PDF_minus[::-1], label = r'$\tilde{J}_{+}$', linestyle = 'dotted', color = "red", zorder = 1)


# axs[1,1].legend(loc = 'lower right')



axs[0,0].set_ylabel(r'$U/U_0$')

axs[0,1].yaxis.tick_right()
axs[0,1].yaxis.set_label_position("right")


# axs[1,0].set_xlim((-0.01, 0.01))
axs[1,0].set_xlabel(r'$\varphi/L$')
axs[1,0].set_ylabel(r'$\widetilde{J}_{\pm}$')
axs[1,0].set_yscale('log')


axs[1,1].set_xlabel(r'$\varphi /L$')
# axs[1,1].set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
axs[1,1].set_yscale('log')
axs[1,1].yaxis.tick_right()
axs[1,1].yaxis.set_label_position("right")


axs[0,0].set_yticks(
    ticks=[-4, 0, 4],
    labels=[r'$-4$', r'$0$', r'$4$']
)

axs[0,1].set_yticks(
    ticks=[-4, 0, 4],
    labels=[r'$-4$', r'$0$', r'$4$']
)


axs[1,0].set_xticks(
    ticks=[-0.5, 0, 0.5],
    labels=[r'$-\pi$', r'$0$', r'$\pi$']
)

axs[1,1].set_xticks(
    ticks=[-0.5, 0, 0.5],
    labels=[r'$-\pi$', r'$0$', r'$\pi$']
)


axs[1,0].set_yticks(
    ticks=[10**-2, 10**0, 10**2],
    labels=[r'$10^{-2}$', r'$10^{0}$', r'$10^2$']
)

axs[1,1].set_yticks(
    ticks=[10**-4, 10**-2],
    labels=[r'$10^{-4}$', r'$10^{-2}$']
)




import matplotlib.patches as patches
filled_circle = patches.Ellipse((0.47, 0.71), width=0.04, height=0.04, color='black', transform=axs[0,0].transAxes)
axs[0,0].add_patch(filled_circle)

empty_circle = patches.Ellipse((0.53, 0.32),width=0.04, height=0.04, color='red', transform=axs[0,0].transAxes)
axs[0,0].add_patch(empty_circle)



'''
axs[0,0].text(0.3, 0.78, r'$U_{-}$', transform=axs[0, 0].transAxes,
         fontsize=11, va='top', ha="right")
axs[0,0].text(0.3, 0.29, r'$U_{+}$', transform=axs[0, 0].transAxes,
         fontsize=11, va='top', ha="right")

axs[0,1].text(0.3, 0.83, r'$U_{-}$', transform=axs[0, 1].transAxes,
         fontsize=11, va='top', ha="right")
axs[0,1].text(0.3, 0.25, r'$U_{+}$', transform=axs[0, 1].transAxes,
         fontsize=11, va='top', ha="right")

axs[1,0].text(0.38, 0.3, r'$\tilde{J}_{-}$', transform=axs[1, 0].transAxes,
         fontsize=11, va='top', ha="right")
axs[1,0].text(0.73, 0.3, r'$\tilde{J}_{+}$', transform=axs[1, 0].transAxes,
         fontsize=11, va='top', ha="right")

axs[1,1].text(0.46, 0.9, r'$\tilde{J}_{-}$', transform=axs[1, 1].transAxes,
         fontsize=11, va='top', ha="right")
axs[1,1].text(0.63, 0.9, r'$\tilde{J}_{+}$', transform=axs[1, 1].transAxes,
         fontsize=11, va='top', ha="right")
'''

plt.tight_layout()


plt.savefig('TUR_Ratchet_a_0_L_PDF.svg')
plt.show()



