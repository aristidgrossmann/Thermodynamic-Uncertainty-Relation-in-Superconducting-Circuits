import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
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



############################     LR CIRCUIT     #########################################

def analytical_Voltage_mean(i):
    return np.where(i <= 1, 0, np.sqrt(i**2 -1))


i_0 = np.linspace(0, 3.1, 1000)
fig, ax = plt.subplots(figsize = (4.5, 3))

ax.plot(i_0, analytical_Voltage_mean(i_0), linestyle = '-', color = 'black')
ax.plot(i_0, i_0, linestyle = '--', color = 'black', alpha = 0.5)
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

plt.tight_layout()
plt.savefig('JJ_no_noise_I_V.pdf')
plt.show()


