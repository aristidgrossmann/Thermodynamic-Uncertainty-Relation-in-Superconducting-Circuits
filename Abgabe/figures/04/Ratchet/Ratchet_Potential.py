import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/Ratchet')

sys.path.append(os.getcwd())

import Exact_Solutions 

Exact_Solutions.L = 1

Exact_Solutions.V_0 = 1



Exact_Solutions.V = Exact_Solutions.V_Ratchet

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

# Define phi range (0 to 2π)
phi = np.linspace(0.1, 2.9, 1000)

# Define i values
i_values = [2,4,6]  #6,4,2

a_vals = [1, 0.5, 1e-4]
a_labels = ["L", "L/2", "0"]

# Plot
fig, ax = plt.subplots(figsize=(4.5, 3))
colors = plt.cm.viridis(np.linspace(0.0,0.7,len(i_values)))
colors = ["black", "blue", "red"]

for i in range(len(i_values)):
    V_vals = np.ones(len(phi))
    Exact_Solutions.a = a_vals[i]* Exact_Solutions.L
    Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
    for j in range(len(phi)):
        V_vals[j] = Exact_Solutions.V(phi[j], -i_values[i], 1)
    a_vals[0] = 0
    ax.plot(phi, V_vals, label=f"$i_0/U_0 = {i_values[i]}, a = {a_labels[i]}$", color = colors[i])

# Customize plot
ax.set_xlabel(r'$\varphi/L$')
ax.set_ylabel(r'$U/U_0$')
ax.set_xticks(
    ticks=[0, 1, 2, 3],
    labels=[r'$0$', r'$1$', r'$2$', r'$3$']
)
ax.set_yticks(
    ticks=[0, -5, -10, -15],
    labels=[r'$0$', r'$-5$', r'$-10$', r'$-15$']
)
plt.legend()
plt.tight_layout()
plt.savefig('Sawtooth_Potential.svg')
plt.show()