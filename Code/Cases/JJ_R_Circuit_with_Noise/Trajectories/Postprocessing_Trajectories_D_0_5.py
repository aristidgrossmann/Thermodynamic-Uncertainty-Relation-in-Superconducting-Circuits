import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.integrate import quad


def MFPT(phi, i_0, D):

    L = 2*np.pi
    mean_velocity = np.zeros_like(phi)

    if D <= 1e-3:
        D = 1e-3

    for i, phi_val in enumerate(phi):

        def integrand(alpha):
            return np.exp((- i_0*alpha + np.cos(phi_val-alpha) - np.cos(phi_val))/D)
        
        integral, _ = quad(integrand, 0, L)
        denominator = D*(1-np.exp(-L*i_0/D))
        mean_velocity[i] =  denominator/integral

    return mean_velocity
    
def T_period_no_noise(i_0):
    return 2*np.pi/np.sqrt(i_0**2 -1)

i_0 = [1.05, 1.36667, 1.6833333, 2]
D = [0.5, 0.5, 0.5, 0.5]
colors = plt.cm.viridis(np.linspace(0,0.9,len(i_0)))

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_with_Noise/Trajectories')


time_0 = np.genfromtxt('0__time.txt')
mean_0 = np.genfromtxt('0__mean.txt')
variance_0 = np.genfromtxt('0__variance.txt')
mean_velocity_0 = np.genfromtxt('0__mean_velocity.txt')

time_1 = np.genfromtxt('1__time.txt')
mean_1 = np.genfromtxt('1__mean.txt')
variance_1 = np.genfromtxt('1__variance.txt')
mean_velocity_1 = np.genfromtxt('1__mean_velocity.txt')

time_2 = np.genfromtxt('2__time.txt')
mean_2 = np.genfromtxt('2__mean.txt')
variance_2 = np.genfromtxt('2__variance.txt')
mean_velocity_2 = np.genfromtxt('2__mean_velocity.txt')

time_3 = np.genfromtxt('3__time.txt')
mean_3 = np.genfromtxt('3__mean.txt')
variance_3 = np.genfromtxt('3__variance.txt')
mean_velocity_3 = np.genfromtxt('3__mean_velocity.txt')


# relative_current_no_Noise = np.linspace(np.min(Relative_Current), np.max(Relative_Current), 1000)


fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1 ,figsize = (6, 9))

ax1.plot(time_0/T_period_no_noise(i_0[0]), mean_0/(2*np.pi), color = colors[0])
ax1.plot(time_1/T_period_no_noise(i_0[1]), mean_1/(2*np.pi), color = colors[1])
ax1.plot(time_2/T_period_no_noise(i_0[2]), mean_2/(2*np.pi), color = colors[2])
ax1.plot(time_3/T_period_no_noise(i_0[3]), mean_3/(2*np.pi), color = colors[3])



ax1.grid(True)
ax1.set_xlabel('time')
ax1.set_ylabel('mean')



ax2.plot(time_0/T_period_no_noise(i_0[0]), variance_0, color = colors[0])
ax2.plot(time_1/T_period_no_noise(i_0[1]), variance_1, color = colors[1])
ax2.plot(time_2/T_period_no_noise(i_0[2]), variance_2, color = colors[2])
ax2.plot(time_3/T_period_no_noise(i_0[3]), variance_3, color = colors[3])


ax2.grid(True)
ax2.set_xlabel('time')
ax2.set_ylabel('variance')


ax3.plot(time_0/T_period_no_noise(i_0[0]), variance_0/mean_0 *i_0[0]/D[0], color = colors[0])
ax3.plot(time_1/T_period_no_noise(i_0[1]), variance_1/mean_1 *i_0[1]/D[1], color = colors[0])


# ax3.plot(mean_0, mean_velocity_0, color = colors[0])
# ax3.plot(mean_0, MFPT(mean_0, i_0[0], D[0]), linestyle = '--', color = colors[0])

# ax3.plot(mean_1, mean_velocity_1, color = colors[1])


ax3.grid(True)
ax3.set_xlabel('time')
ax3.set_ylabel('TUR')

plt.tight_layout()
plt.show()

# ax.plot(Relative_Current, Voltage_mean, marker = 'x', s = 5, color = 'red', label = 'numerical, with noise')
# ax.fill_between(Relative_Current, Voltage_mean - Voltage_mean_std, Voltage_mean + Voltage_mean_std, alpha = 0.3, color = 'red', edgecolor = 'none')
# ax.plot(relative_current_no_Noise, analytical_Voltage_mean_no_Noise(relative_current_no_Noise), color = 'blue', label = 'analytical, no noise')
# ax.plot(Relative_Current, V_vals_stable_2, label = 'analytical, with noise', color = 'black')
# ax.set_title(r'I-V curve for $\sqrt{\frac{2 k_{\mathrm{B}}T}{E_{\mathrm{J}}}} = 10$')
# ax.set_xlabel(r"$I_0 / I_{\mathrm{c}}$", fontsize=12)
# ax.set_ylabel(r"$\overline{\langle V \rangle } \, 2e/\hbar$", fontsize=12)
# ax.grid(True)
# ax.legend()


# plt.savefig('JJ_R_Circuit_Noise_IV_ratio_10.pdf')
# plt.show()


