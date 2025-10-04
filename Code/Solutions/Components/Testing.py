import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/Components')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, Unc_prod, I_tilde_plus_minus, I_tilde_minus_plus



def central_differences(x, y):
  
    dy_dx = np.zeros_like(y)
    
    # Central differences for interior points
    dy_dx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    
    # Forward difference at the start
    dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])
    
    # Backward difference at the end
    dy_dx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    
    return dy_dx


i_0_min = 0.0
i_0_max = 2
N_i_0 = 58
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)


D_vals = np.array([0.5, 1, 2, 5, 10])
D_labels = [r"$D = 0.5$", r"$D = 1$", r"$D = 2$", r"$D = 5$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))



fig, ax = plt.subplots()

for i in range(len(D_vals)):
    Test = np.ones_like(i_0_vals)
    for j in range(len(i_0_vals)):
        Test[j] = Unc_prod(i_0 = i_0_vals[j], D = D_vals[i]) 


    ax.plot(i_0_vals, Test, label=D_labels[i], color = colors[i+1])





ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
plt.tight_layout()
plt.legend()
# plt.savefig('Test.pdf')
plt.show()



