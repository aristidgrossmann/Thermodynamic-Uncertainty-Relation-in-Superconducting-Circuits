import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import  L, Probability_Current_double_integral




i_0_min = 1.25
i_0_max = 2
N_i_0 = 51
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.5, 1e0, 10])
D_labels = [r"$D = 0.5$", r"$D = 1$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))


fig, ax = plt.subplots()

for i in range(len(D_vals)):
    V_vals = np.ones_like(i_0_vals)
    for j in range(len(V_vals)):
        V_vals[j] = L*Probability_Current_double_integral(i_0 = i_0_vals[j], D = D_vals[i])

    ax.plot(i_0_vals, V_vals, label=D_labels[i], color = colors[i+1])




ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Mean_Velocity.pdf')
plt.show()



