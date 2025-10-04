import numpy as np
from matplotlib import pyplot as plt
import math
import os
from scipy.integrate import quad
from Exact_Solutions import L, analytical_Voltage_mean_no_Noise, Probability_Current_double_integral , Delta_T_2_paper, T_1_paper

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Solutions')

#########################   Helper Functions   ###################################
def omega_0_calc(i_0):
     return np.sqrt(i_0**2 -1.0)

def theta_2_dot_time_dependent(tau, i_0, D):
        omega_0 = omega_0_calc(i_0)
        pre_exp = (2 * D / omega_0**2) * np.cos(omega_0 * tau) * (
            np.sin(omega_0 * tau) + i_0
        )
        return pre_exp * np.exp(-0.5*variance_theta_time_dependent(tau, i_0, D))

def variance_theta_time_dependent(tau, i_0, D):
        omega_0 = omega_0_calc(i_0)
        return -D/omega_0**2*((np.sin(2*omega_0*tau)/4.0 +2.0*i_0*np.cos(omega_0*tau)-2.0*i_0)/omega_0
                             -(i_0**2 +1.0/2.0)*tau)

#########################   Mean Velocity   ###################################
def theta_2_dot_time_average(i_0, D):

    omega_0 = omega_0_calc(i_0)

    T = 2 * np.pi / omega_0
    average, _ = quad(theta_2_dot_time_dependent, 0, T, args = (i_0, D))
    average /= T
    return average

def mean_velocity_approximate(i_0, D):   
    return np.sqrt(i_0**2 -1) + theta_2_dot_time_average(i_0, D)


#########################   Variance velocity   ###################################
def variance_time_dependent(tau, i_0, D):
     omega_0 = omega_0_calc(i_0)
     prefactor = 1/(np.sin(omega_0*tau) + i_0)**2
     return prefactor *variance_theta_time_dependent(tau, i_0, D)


def variance_velocity_approximate(i_0, D):
    omega_0 = omega_0_calc(i_0)
    T = 2 * np.pi / omega_0

    average, _ = quad(variance_time_dependent, 0, T, args = (i_0, D))
    average /= T
    return average

    



i_0_min = 1.1
i_0_max = 2
N_i_0 = 10

i_0_vals = np.linspace(i_0_min, i_0_max, N_i_0)
D_vals = np.array([1e-2, 1e-1])
D_labels = [r"$D = 1\cdot 10^{-2}$", r"$D = 1\cdot 10^{-1}$"]
colors = plt.cm.viridis(np.linspace(0,0.8,len(D_vals)+1))


fig, ax = plt.subplots(figsize = (6, 4.5))
ax.plot(i_0_vals, analytical_Voltage_mean_no_Noise(i_0_vals), label = r'D = 0', color = colors[0])

for i in range(1,len(D_vals)):
    V_vals_analytical = np.ones_like(i_0_vals)
    V_vals_approx = np.ones_like(i_0_vals)
    for j in range(len(V_vals_analytical)):
        V_vals_analytical[j] = L*Probability_Current_double_integral(i_0 = i_0_vals[j], D = D_vals[i])
        V_vals_approx[j] = mean_velocity_approximate(i_0_vals[j], D_vals[i])

    ax.plot(i_0_vals, V_vals_analytical, label=D_labels[i], color = colors[i+1])
    ax.plot(i_0_vals, V_vals_approx, label=D_labels[i] + ' approx', color = 'red')

ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
ax.grid(True)
plt.legend()
plt.savefig('Mod_Time_I-V.pdf')
plt.show()



#############################   Variance   ############################
fig, ax = plt.subplots(figsize = (6, 4.5))

for i in range(1,len(D_vals)):
    V_ariance_vals_approx = np.ones_like(i_0_vals)
    Variance_vals_analytical = np.ones_like(i_0_vals)
    for j in range(len(Variance_vals_analytical)):
        Variance_vals_analytical[j] = L**2*Delta_T_2_paper(i_0 = i_0_vals[j], D = D_vals[i]) / T_1_paper(i_0 = i_0_vals[j], D = D_vals[i])**3
        V_ariance_vals_approx[j] = variance_velocity_approximate(i_0 = i_0_vals[j], D = D_vals[i])

    ax.plot(i_0_vals, Variance_vals_analytical, label=D_labels[i], color = colors[i+1])
    ax.plot(i_0_vals, V_ariance_vals_approx, label=D_labels[i] + ' approx', color = 'red')

ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
ax.set_ylabel("Variance velo", fontsize=12)
ax.grid(True)
plt.legend()
plt.savefig('Variance_Velo_approx.pdf')
plt.show()



