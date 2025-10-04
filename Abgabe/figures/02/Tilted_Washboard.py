import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/02')

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
phi = np.linspace(0, 4*np.pi, 1000)

# Define i values
i_values = [0.2, 1, 2]

# Plot
fig, ax = plt.subplots(figsize=(4.5, 3))
colors = plt.cm.viridis(np.linspace(0.0,0.7,len(i_values)))
colors = ["black", "blue", "red"]


for i in range(len(i_values)):
    ax.plot(phi, -i_values[i]*phi - np.cos(phi), label=f"$i_0 = {i_values[i]}$", color = colors[i])

# Customize plot
ax.set_xlabel(r'$\varphi$')
ax.set_ylabel(r'$U$')
ax.set_xticks(
    ticks=[0, 2*np.pi, 4*np.pi],
    labels=[r'$0$', r'$2\pi$', r'$4\pi$']
)
ax.set_yticks(
    ticks=[0, -2*np.pi, -4*np.pi, -6*np.pi, -8*np.pi],
    labels=[r'$0$', r'$-2\pi$', r'$-4\pi$', r'$-6\pi$', r'$-8\pi$']
)





plt.legend()
plt.tight_layout()
plt.savefig('Tilted_Washboard_Potential.svg')
plt.show()