import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os
import sys


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/Vortrag/figs/02')
sys.path.append( os.getcwd())

from Exact_Solutions import Mean_Velocity




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


colors = ["black", "blue", "red"]


def analytical_Voltage_mean(i):
    return np.where(i <= 1, 0, np.sqrt(i**2 -1))


i_0_min = 0.01
i_0_max = 3.0
N_i_0 = 1000
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.1, 0.5, 1.5])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1.5$"]

fig, ax = plt.subplots(figsize = (4.5, 3))

ax.plot(i_0_vals, analytical_Voltage_mean(i_0_vals), linewidth = 1.5, linestyle = 'dashed', color = 'purple', zorder = 1)
ax.plot(i_0_vals, i_0_vals, linestyle = 'dashed', linewidth = 2, color = 'black', zorder = 1)

for i in range(len(D_vals)):
    V_vals = np.ones_like(i_0_vals)

    for j in range(len(V_vals)): 
        V_vals[j] = Mean_Velocity(i_0 = i_0_vals[j], D = D_vals[i])
        print(j)

    # ax.plot(i_0_vals, Uncertainty_product_vals, label=D_labels[i], linestyle = '-', color = colors[i], zorder = 1)
    ax.plot(i_0_vals, V_vals, linestyle = '-', label = D_labels[i], color = colors[i], zorder = 2)


ax.set_xlabel('$I_0 / I_{\mathrm{c}}$')
ax.set_ylabel(r'$\overline{V}/(R I_{\mathrm{c}})$')
ax.set_xticks([0,1,2,3])
ax.set_yticks([0,1,2,3])

slope = 0.6
x_text = 1.5
y_text = slope* x_text

# Compute the angle of the line (in degrees)
angle_rad = np.arctan(slope)
angle_deg = np.degrees(angle_rad)

# Add text above the line, offset slightly in y
ax.text(x_text, y_text + 0.6, "$\overline{V} = R I_0$", rotation=angle_deg,
        rotation_mode='anchor', ha='center', va='bottom')

ax.legend(loc = 'lower right')

plt.tight_layout()
plt.savefig('JJ_noise_I_V.svg')
plt.show()


