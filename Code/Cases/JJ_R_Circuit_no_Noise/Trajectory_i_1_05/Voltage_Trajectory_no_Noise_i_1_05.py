import numpy as np
from matplotlib import pyplot as plt
import os


os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Cases/JJ_R_Circuit_no_Noise/Trajectory_i_1_05')

#read files: 
time = np.genfromtxt("time.txt")
MC_mean = np.genfromtxt("mean.txt")

i = 1.05
T = 2*np.pi/np.sqrt(i**2 -1)

def central_difference_derivative(f, t):
    """
    Compute the numerical derivative of f with respect to t using central differences.

    Parameters:
    f (numpy.ndarray): Array of function values.
    t (numpy.ndarray): Array of corresponding independent variable values (e.g., time or position).

    Returns:
    numpy.ndarray: The derivative df/dt computed via central differences.
                  The first and last points use forward and backward differences, respectively.
    """
    if len(f) != len(t):
        raise ValueError("Arrays f and t must have the same length.")
    if len(f) < 2:
        raise ValueError("Arrays must have at least 2 elements to compute a derivative.")

    dfdt = np.ones_like(f)
    dt = np.diff(t)  # Differences between consecutive t values

    # Central differences for interior points
    dfdt[1:-1] = (f[2:] - f[:-2]) / (t[2:] - t[:-2])

    # Forward difference for the first point
    dfdt[0] = (f[1] - f[0]) / dt[0]

    # Backward difference for the last point
    dfdt[-1] = (f[-1] - f[-2]) / dt[-1]

    return dfdt

def Voltage_Analytical(t_values, i):
    """
    Compute the expression for an array of t values and a given i.
    
    Parameters:
    t_values (array-like): Array of input values for t.
    i (float): Positive real number (i > 1 to ensure real-valued results).
    
    Returns:
    numpy.ndarray: Array of computed values for the given expression.
    """
    t = np.asarray(t_values, dtype=float)
    sqrt_i2_minus_1 = np.sqrt(i**2 - 1)
    
    numerator = i * (i**2 - 1)
    
    denominator_part1 = i * (np.sqrt(1 - 1/i**2) * np.sin(sqrt_i2_minus_1 * t) + i)
    denominator_part2 = np.cos(sqrt_i2_minus_1 * t)
    denominator = denominator_part1 + denominator_part2
    
    # Avoid division by zero (though unlikely for i > 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = numerator / denominator
    
    return result



############################     LR CIRCUIT     #########################################

def analytical_Voltage_mean(i):
    return np.where(i <= 1, 0, np.sqrt(i**2 -1))


fig, ax = plt.subplots(figsize = (6, 4.5))

ax.plot(time/T, central_difference_derivative(MC_mean, time)/(2*np.pi), color = 'red', label = 'numerical')
ax.plot(time/T, Voltage_Analytical(time - T/2, i)/(2*np.pi), color = 'blue', label = 'analytical')
#ax.plot(Relative_Current, analytical_Voltage_mean(Relative_Current), color = 'blue', label = 'analytical')
ax.set_title(r'Trajectory for $i = 1.05$')
ax.set_xlabel(r'$t /(t_0 \hat{T})$')
ax.set_ylabel(r'$V t_0/(2\pi \phi_0) $ ')
ax.grid(True)
ax.legend(loc = 'upper right')



plt.tight_layout()
plt.savefig('JJ_R_Circuit_Voltage_Trajectory_i_1_05.pdf')
plt.show()


