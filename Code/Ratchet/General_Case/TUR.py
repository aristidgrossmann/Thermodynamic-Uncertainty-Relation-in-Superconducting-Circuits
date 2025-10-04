import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Ratchet/General_Case')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

import Exact_Solutions 

Exact_Solutions.L = 1
Exact_Solutions.a = 0.3* Exact_Solutions.L
Exact_Solutions.b = Exact_Solutions.L - Exact_Solutions.a
Exact_Solutions.V_0 = 2


Exact_Solutions.V = Exact_Solutions.V_Ratchet


i_0_min = 1
i_0_max = 10
N_i_0 = 40
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([0.1, 0.5, 1, 2, 5, 10])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1$", r"$D = 2$", r"$D = 5$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))

# plt.figure()
# x = np.linspace(-2, 2, 100)
# y = np.ones_like(x)
# for i in range(len(x)):
#     y[i] = Exact_Solutions.V_Ratchet(x[i], 1, 0.2)
# plt.figure()
# plt.plot(x, y)
# plt.show()


fig, ax = plt.subplots(figsize = (7, 5.2))

for i in range(0,len(D_vals)):
    TUR_vals = np.ones_like(i_0_vals)
    TUR_vals_explicit_ratchet = np.ones_like(i_0_vals)
    for j in range(len(TUR_vals)): 
        TUR_vals[j] = Exact_Solutions.Uncertainty_Product_Ratchet(i_0 = i_0_vals[j], D = D_vals[i])
        # TUR_vals_explicit_ratchet[j] = Exact_Solutions.T_1_Ratchet(i_0 = i_0_vals[j], D = D_vals[i])

    ax.plot(i_0_vals, TUR_vals, label=D_labels[i], color = colors[i+1])
    # ax.plot(i_0_vals, TUR_vals_explicit_ratchet,label=D_labels[i], color = colors[i+1])




ax.hlines(2, i_0_min - 0.05, i_0_max+0.05, linestyle = '--', linewidth = 2, label = 'TUR lower bound', color = 'black')



ax.set_xlabel(r"$i_0$", fontsize=12)
ax.set_ylabel("Uncertainty Product", fontsize=12)
ax.set_yscale('log')

plt.tight_layout()
plt.legend()


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Ratchet/General_Case')

plt.savefig('Ratchet_General_Case_TUR_analytical.pdf')
plt.show()



