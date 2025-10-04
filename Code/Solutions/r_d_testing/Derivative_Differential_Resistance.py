import numpy as np
from matplotlib import pyplot as plt
import os
import sys

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/r_d_testing')

parent_dir = os.path.abspath("..")
sys.path.append( parent_dir)

from Exact_Solutions import L, r_d, I_tilde_plus_minus, derivative_r_d, Uncertainty_product, T_1_paper, Variance_Velocity



def central_differences(x, y):
  
    dy_dx = np.zeros_like(y)
    
    # Central differences for interior points
    dy_dx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    
    # Forward difference at the start
    dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])
    
    # Backward difference at the end
    dy_dx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    
    return dy_dx


i_0_min = 0.1
i_0_max = 2
N_i_0 = 3
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

N_i_0_numerical = 58
i_0_vals_numerical = np.linspace(i_0_min, i_0_max, N_i_0_numerical)

D_vals = np.array([0.1, 0.5, 1, 2, 5, 10])
D_labels = [r"$D = 0.1$", r"$D = 0.5$", r"$D = 1$", r"$D = 2$", r"$D = 5$", r"$D = 10$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))



fig, ax = plt.subplots(figsize = (7,7))

for i in range(3,len(D_vals)):
    r_d_analytical = np.ones_like(i_0_vals_numerical)
    r_d_derivative = np.ones_like(i_0_vals)
    for j in range(len(i_0_vals)):
        print(j + i*j)
        r_d_derivative[j] = derivative_r_d(i_0 = i_0_vals[j], D = D_vals[i])

    for j in range(len(i_0_vals_numerical)):
        r_d_analytical[j] =  r_d(i_0 = i_0_vals_numerical[j], D = D_vals[i])

    ax.plot(i_0_vals_numerical, central_differences(i_0_vals_numerical, r_d_analytical), label=D_labels[i] + ' r_d central differences', color = colors[i+1])
    ax.scatter(i_0_vals, r_d_derivative, label=D_labels[i] + ' r_d derivative', color = 'red')





ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
plt.tight_layout()
plt.legend()
# plt.savefig('Derivative Differential Resistance.pdf')
plt.show()



