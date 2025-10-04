import numpy as np
import matplotlib.pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/Tilted_Washboard_Potential')

# Define phi range (0 to 2π)
phi = np.linspace(0, 2, 1000)

# Define i values
i_values = [0.2, 1, 2]

# Plot
fig, ax = plt.subplots(figsize=(6, 4.5))
colors = plt.cm.viridis(np.linspace(0.0,0.7,len(i_values)))
for i in range(len(i_values)):
    ax.plot(phi, -i_values[i]*phi - np.cos(phi*2*np.pi)/(2*np.pi), label=f"$i_0$ = {i_values[i]}", color = colors[i])

# Customize plot
ax.set_xlabel(r"$\varphi /(2\pi )$", fontsize=12)
ax.set_ylabel(r"$V/(2\pi)$", fontsize=12)
ax.set_xticks([0,1,2])
ax.set_yticks([0,-1,-2, -3, -4])
plt.legend()
plt.savefig('Tilted_Washboard_Potential.pdf')
plt.show()