import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.integrate import dblquad


def compute_V_stable_2(i_0, D):
    """
    Compute the expression J numerically.
    
    Parameters:
    i_0: Current ratio I_0/I_c
    D : float
        Diffusion coefficient D
    L : float
        Length of the interval [0, L]
    """

    L = 2*np.pi
    # Helper function for the exponent term

    def integrand(x, y):
        return np.exp((np.cos(x)-np.cos(x+y) - i_0*y)/D)
    
    integral, _ = dblquad(
    integrand,
    0, L,  # x limits
    0, L #y limits
    )

    denominator = 2*np.pi*D*(1-np.exp(-L*i_0/D))
    
    J = denominator/integral
    
    return J


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/IV_Curve_Noise_ratio_0_1')

B = 0.1
D = B**2/2

#read files: 
Relative_Current = np.genfromtxt("Relative_Current_IV.txt")

Voltage_mean = np.genfromtxt("Voltage_mean_IV.txt")
Voltage_mean_std = np.genfromtxt("Voltage_mean_std_IV.txt")

def analytical_Voltage_mean_no_Noise(i):
    vals = np.zeros_like(i)
    for j in range(len(i)): 
        if i[j] <= 1:
            continue
        else: 
            vals[j] = np.sqrt(i[j]**2 -1)
    #return np.where(i <= 1, 0, np.sqrt(i**2 -1))
    return vals

relative_current_no_Noise = np.linspace(np.min(Relative_Current), np.max(Relative_Current), 1000)

V_vals_stable_2 = np.zeros_like(Relative_Current)
for i in range(len(Relative_Current)):
    V_vals_stable_2[i] = compute_V_stable_2(Relative_Current[i], D)


fig, ax = plt.subplots(figsize = (6, 4.5))

ax.scatter(Relative_Current, Voltage_mean, marker = 'x', s= 5, color = 'red', label = 'numerical, with noise')
ax.plot(Relative_Current, V_vals_stable_2, label = 'analytical, with noise', color = 'black')
#ax.fill_between(Relative_Current, Voltage_mean - Voltage_mean_std, Voltage_mean + Voltage_mean_std, alpha = 0.3, color = 'red', edgecolor = 'none')
ax.plot(relative_current_no_Noise, analytical_Voltage_mean_no_Noise(relative_current_no_Noise), color = 'blue', label = 'analytical, no noise')

ax.set_title(r'I-V curve for $\sqrt{\frac{2 k_{\mathrm{B}}T}{E_{\mathrm{J}}}} = 0.1$')
ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
ax.legend()



plt.tight_layout()
plt.savefig('JJ_R_Circuit_Noise_IV_ratio_0_1.pdf')
plt.show()


