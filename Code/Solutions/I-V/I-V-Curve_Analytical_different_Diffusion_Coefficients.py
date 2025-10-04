import numpy as np
from matplotlib import pyplot as plt
import math
import os
from Exact_Solutions import analytical_Voltage_mean_no_Noise, analytical_Voltage_Probability_Current_double_integral 


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions/I-V')
i_0_min = 0
i_0_max = 2
N_i_0 = 201
i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)

D_vals = np.array([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e5])
D_labels = [r"$D = 1\cdot 10^{-3}$", r"$D = 1\cdot 10^{-2}$", r"$D = 1\cdot 10^{-1}$", 
            r"$D = 1\cdot 10^{0}$", r"$D = 1\cdot 10^{1}$", r"$D = 1\cdot 10^{5}$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))


fig, ax = plt.subplots()
ax.plot(i_0_vals, analytical_Voltage_mean_no_Noise(i_0_vals), label = r'D = 0', color = colors[0])

for i in range(len(D_vals)):
    V_vals = np.ones_like(i_0_vals)
    for j in range(len(V_vals)):
        V_vals[j] = analytical_Voltage_Probability_Current_double_integral(i_0 = i_0_vals[j], D = D_vals[i])

    ax.plot(i_0_vals, V_vals, label=D_labels[i], color = colors[i+1])

ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
plt.legend()
plt.savefig('I_V_Curve_Analytical_Different_Diffusion_Coefficients.pdf')
plt.show()
