import numpy as np
import math
from scipy.integrate import quad, dblquad
from scipy.special import i0
from matplotlib import pyplot as plt
import os

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions')


infinity = np.inf
infinity = 1e2

L = 2*np.pi

def V(phi, i_0, D):
    return (-i_0*phi - np.cos(phi))/D

def I_tilde_plus_minus(x, i_0, D):

    def integrand(alpha, i_0, D):
        return np.exp(- V(alpha, i_0, D))
    
    integral, _ = quad(integrand, x - L, x, args= (i_0, D))
    return integral*np.exp( V(x, i_0, D))

def I_tilde_minus_plus(x, i_0, D):

    def integrand(alpha, i_0, D):
        return np.exp( V(alpha, i_0, D))
    
    integral, _ = quad(integrand, x , x + L, args= (i_0, D))
    return integral*np.exp( -V(x, i_0, D))


#####  FPT Mean  ########
def T_1(i_0, D):
    result, _ = quad(I_tilde_plus_minus, 0, L, args= (i_0, D))
    return result/(D*(1 - np.exp(- i_0*L/D)))

#####  FPT Variance  ########
def Delta_T_2(i_0, D):
    L = 2*np.pi
    
    def integrand_z(z, i_0, D):
        return I_tilde_plus_minus(z, i_0, D)**2 * I_tilde_minus_plus(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return 2/D**2 /(1 - np.exp(- i_0*L/D)**3) * result  


def central_differences(x, y):
  
    dy_dx = np.zeros_like(y)
    
    # Central differences for interior points
    dy_dx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    
    # Forward difference at the start
    dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0])
    
    # Backward difference at the end
    dy_dx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    
    return dy_dx



######   I-V curve with MFPT #####################################################
i_0_min = 1
i_0_max = 2
N_i_0 = 21


i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)
D_vals = np.array([1e-1, 1e0, 1e1, 1e5])
D_labels = [r"$D = 1\cdot 10^{-1}$", r"$D = 1\cdot 10^{0}$", r"$D = 1\cdot 10^{1}$", r"$D = 1\cdot 10^{5}$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))


fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1 ,figsize = (6, 9))

for i in range(len(D_vals)):
    V_vals = np.ones_like(i_0_vals)
    TUR_vals = np.ones_like(i_0_vals)
    for j in range(len(V_vals)):
        V_vals[j] = 2*np.pi/T_1(i_0 = i_0_vals[j], D = D_vals[i])
        TUR_vals[j] = 2*np.pi*Delta_T_2(i_0 = i_0_vals[j], D = D_vals[i])/T_1(i_0 = i_0_vals[j], D = D_vals[i])**2 *i_0_vals[j]/D_vals[i]

    ax1.plot(i_0_vals, V_vals, label=D_labels[i], color = colors[i+1])
    ax2.plot(i_0_vals, TUR_vals, label=D_labels[i], color = colors[i+1])
    ax3.plot(i_0_vals, central_differences(i_0_vals, V_vals)*i_0_vals*2*np.pi/V_vals, label=D_labels[i], color = colors[i+1])





ax1.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax1.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax1.grid(True)
ax1.legend()

ax2.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax2.set_ylabel('TUR LHS', fontsize=12)
ax2.grid(True)
ax2.legend()

ax3.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax3.set_ylabel(r"$r_d$", fontsize=12)
ax3.grid(True)
ax3.legend()

plt.tight_layout()
plt.savefig('TUR.pdf')
plt.show()