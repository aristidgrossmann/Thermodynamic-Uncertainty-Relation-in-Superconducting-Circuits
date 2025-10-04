import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/r_d_testing')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, r_d, I_tilde_plus_minus, Uncertainty_product, T_1_paper, Mean_Velocity, Variance_Velocity, I_tilde_plus_minus, I_tilde_minus_plus



def central_differences(x, y):
  
    dy_dx = np.zeros_like(y)
    
    # Central differences for interior points
    dy_dx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    
    # Forward difference at the start
    dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])
    
    # Backward difference at the end
    dy_dx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    
    return dy_dx


i_0_min = 0.5
i_0_max = 2
N_i_0 = 4
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D = 0.1

noof_x_vals = 40
x_vals = np.linspace(0, L, noof_x_vals)
colors = plt.cm.viridis(np.linspace(0,0.8,len(i_0_vals)+1))



fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize = (8,10))

for i in range(len(i_0_vals)):
    test_vals = np.ones_like(x_vals)


    for j in range(len(x_vals)):
        test_vals[j] = 1/(D*(1-np.exp(-i_0_vals[i]*L/D)))**2*(I_tilde_minus_plus(x_vals[j], i_0 = i_0_vals[i], D = D)* I_tilde_plus_minus(x_vals[j], i_0 = i_0_vals[i], D = D))

    ax1.plot(2*x_vals/L, test_vals,  label = i_0_vals[i], color = colors[i+1])

    for j in range(len(x_vals)):
        test_vals[j] = 1/(D*(1-np.exp(-i_0_vals[i]*L/D)))*I_tilde_plus_minus(x_vals[j], i_0 = i_0_vals[i], D = D)
    ax2.plot(2*x_vals/L, test_vals,  label = i_0_vals[i], color = colors[i+1])

    for j in range(len(x_vals)):
        test_vals[j] = 1/(D*(1-np.exp(-i_0_vals[i]*L/D)))*I_tilde_minus_plus(x_vals[j], i_0 = i_0_vals[i], D = D)
    ax3.plot(2*x_vals/L, test_vals,  label = i_0_vals[i], color = colors[i+1])








ax1.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax1.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax1.grid(True)

ax2.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax2.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax2.grid(True)
ax2.set_title(r"$I_+$")

ax3.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax3.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax3.set_title(r"$I_-$")
ax3.grid(True)
plt.tight_layout()
plt.legend()
plt.show()



