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

from Exact_Solutions import T_1_Ratchet_a_0, T_1_Ratchet_a_L, I_plus_minus_Ratchet_a_0, I_minus_plus_Ratchet_a_0, I_plus_minus_Ratchet_a_L, I_minus_plus_Ratchet_a_L


import Exact_Solutions 

Exact_Solutions.L = 1
Exact_Solutions.V_0 = 1




i_0 = 6
N_phi = 1000
phi_vals = np.linspace(-0.5, 0.5, N_phi)



D_vals = np.array([0.1, 0.3, 0.5])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]




# fig, (ax1, ax2) = plt.subplots(2, 1,   figsize=(4.5, 5.4), sharex=True)
fig, axs = plt.subplots(2, 2,   figsize=(6.20, 4.5), sharex='col')

### ax11
Exact_Solutions.a = Exact_Solutions.L*0.0
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
Exact_Solutions.V = Exact_Solutions.V_Ratchet
eps = 1e-3
phi_V_vals = np.linspace(-0.5 - eps, 0.5 + eps, 1000)
V_vals = np.ones_like(phi_V_vals)
for j in range(len(V_vals)):
    V_vals[j] = Exact_Solutions.V(phi_V_vals[j], i_0, 1)
axs[0,0].plot(phi_V_vals, V_vals,   color = "purple")
axs[0,0].plot(phi_V_vals, -V_vals -V_vals[0], linestyle = 'dotted',  color = "purple")

### ax12
Exact_Solutions.a = Exact_Solutions.L*1.0
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
Exact_Solutions.V = Exact_Solutions.V_Ratchet
for j in range(len(V_vals)):
    V_vals[j] = Exact_Solutions.V(phi_V_vals[j], i_0, 1)
axs[0,1].plot(phi_V_vals, V_vals,   color = "purple")
axs[0,1].plot(phi_V_vals, -V_vals - V_vals[0], linestyle = 'dotted',  color = "purple")

for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)

    for j in range(len(PDF_minus)): 
        PDF_minus[j] = I_minus_plus_Ratchet_a_0(phi_vals[j], i_0 = i_0, D = D_vals[i])/T_1_Ratchet_a_0(i_0 = i_0, D = D_vals[i])

        if j % 10 == 0:
            print(j)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    axs[1,0].plot(phi_vals, PDF_minus, linestyle = '-', color = colors[i], zorder = 1)
    axs[1,0].plot(phi_vals, PDF_minus[::-1], linestyle = 'dotted', color = colors[i], zorder = 1)

    filename = "PDF_minus_a_0_" + str(D_vals[i]) + ".txt"
    np.savetxt(filename, PDF_minus)




for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)

    for j in range(len(PDF_minus)): 
        PDF_minus[j] = I_minus_plus_Ratchet_a_L(phi_vals[j], i_0 = i_0, D = D_vals[i])/T_1_Ratchet_a_L(i_0 = i_0, D = D_vals[i])
        
        if j % 10 == 0:
            print(j)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    axs[1,1].plot(phi_vals, PDF_minus, linestyle = '-', color = colors[i], zorder = 1)
    axs[1,1].plot(phi_vals, PDF_minus[::-1], linestyle = 'dotted', color = colors[i], zorder = 1)

    filename = "PDF_minus_a_L_" + str(D_vals[i]) + ".txt"
    np.savetxt(filename, PDF_minus)




axs[0,0].set_ylabel(r'$U/U_0$')

axs[0,1].yaxis.tick_right()
axs[0,1].yaxis.set_label_position("right")



axs[1,0].set_xlabel(r'$\varphi/U_0$')
axs[1,0].set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
axs[1,0].set_yscale('log')


axs[1,1].set_xlabel(r'$\varphi / U_0$')
# axs[1,1].set_ylabel(r'$\tilde{\mathcal{Q}}, \tilde{\mathcal{L}}$')
axs[1,1].set_yscale('log')
axs[1,1].yaxis.tick_right()
axs[1,1].yaxis.set_label_position("right")


# ax1.set_xticks(
#     ticks=[10**-2, 10**-1, 10**0, 10**1, 10**2],
#     labels=[r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$']
# )

# ax1.set_yticks(
#     ticks=[2, 10**1, 10**2, 10**3],
#     labels=[r'$2$', r'$10^1$', r'$10^2$', r'$10^3$']
# )


# ax2.set_xticks(
#     ticks=[10**-2, 10**-1, 10**0, 10**1, 10**2],
#     labels=[r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$']
# )

# ax2.set_yticks(
#     ticks=[2, 10],
#     labels=[r'$2$', r'$10$']
# )

# ax1.xaxis.set_minor_locator(LogLocator(subs='all'))  # Auto minor ticks for log scale
# ax1.yaxis.set_minor_locator(LogLocator(subs='all'))
# ax1.tick_params(axis='both', which='minor', labelleft=False, labelbottom=False)  # Hide labels

# ax2.xaxis.set_minor_locator(LogLocator(subs='all'))  # Auto minor ticks for log scale
# ax2.yaxis.set_minor_locator(LogLocator(subs='all'))
# ax2.tick_params(axis='both', which='minor', labelleft=False, labelbottom=False)  # Hide labels

leftorright = "left"
axs[0,0].text(0.98, 0.97, '(a)', transform=axs[0, 0].transAxes,
         fontsize=11, va='top', ha=leftorright)

axs[0,1].text(0.97, 0.97, '(b)', transform=axs[0,1].transAxes,
         fontsize=11, va='top', ha=leftorright)

axs[1,0].text(0.97, 0.97, '(c)', transform=axs[1,0].transAxes,
         fontsize=11, va='top', ha=leftorright)

axs[1,1].text(0.97, 0.97, '(d)', transform=axs[1,1].transAxes,
         fontsize=11, va='top', ha=leftorright)

plt.tight_layout()


plt.savefig('TUR_Ratchet_a_0_PDF.pdf')
plt.show()



