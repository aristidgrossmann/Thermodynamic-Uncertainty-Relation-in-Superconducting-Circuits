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


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/JJ')

sys.path.append( os.getcwd())

from Exact_Solutions import V


i_0 = 0.8
N_phi = 1000
phi_vals = np.linspace(-np.pi/2, 1.5*np.pi, N_phi)



# D_vals = np.array([0.1, 0.3, 0.5])
D_vals = np.array([0.1, 0.5, 1.5])
D_vals = np.array([0.1, 1.5])


D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
colors = plt.cm.viridis(np.linspace(0,0.7,len(D_vals)))

colors = ["black", "blue", "red"]
colors = ["black", "red"]





# fig, (ax1, ax2) = plt.subplots(2, 1,   figsize=(4.5, 5.4), sharex=True)
fig, axs = plt.subplots(2, 1,   figsize=(4.5, 3.5), sharex='col')


phi_V_vals = np.linspace(np.min(phi_vals), np.max(phi_vals), 1000)
V_vals = np.ones_like(phi_V_vals)
for j in range(len(V_vals)):
    V_vals[j] = V(phi_V_vals[j], i_0, 1)
axs[0].plot(phi_V_vals, V_vals,   color = "purple")
axs[0].plot(phi_V_vals, -V_vals +V_vals[-1], linestyle = 'dotted',  color = "purple")


for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)

    filename = "PDF_minus_" + str(D_vals[i]) + ".txt"
    PDF_minus = np.genfromtxt(filename)

    filename = "PDF_plus_" + str(D_vals[i]) + ".txt"
    PDF_plus = np.genfromtxt(filename)

        

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    axs[1].plot(phi_vals, PDF_minus, linestyle = '-', color = colors[i], zorder = 1)
    axs[1].plot(phi_vals, PDF_plus, linestyle = 'dotted', color = colors[i], zorder = 1)



axs[0].set_ylabel(r'$U$')

axs[1].set_xlabel(r'$\varphi$')
axs[1].set_ylabel(r'$\widetilde{J}_{\pm}$')



axs[0].set_yticks(
    ticks=[-4, 0],
    labels=[r'$-4$', r'$0$']
)

axs[1].set_xticks(
    ticks=[0, np.pi],
    labels=[r'$0$', r'$\pi$']
)

axs[1].set_yticks(
    ticks=[0, 1],
    labels=[r'$0$', r'$1$']
)

import matplotlib.patches as patches
filled_circle = patches.Ellipse((0.38, 0.63), width=0.03, height=0.07, color='black', transform=axs[0].transAxes)
axs[0].add_patch(filled_circle)

empty_circle = patches.Ellipse((0.6, 0.44),width=0.03, height=0.07, color='grey', transform=axs[0].transAxes)
axs[0].add_patch(empty_circle)







leftorright = "right"
axs[0].text(0.97, 0.97, '(a)', transform=axs[0].transAxes,
         fontsize=11, va='top', ha=leftorright)

axs[1].text(0.97, 0.97, '(b)', transform=axs[1].transAxes,
         fontsize=11, va='top', ha=leftorright)


plt.tight_layout()


plt.savefig('TUR_JJ_PDF.pdf')
plt.show()



