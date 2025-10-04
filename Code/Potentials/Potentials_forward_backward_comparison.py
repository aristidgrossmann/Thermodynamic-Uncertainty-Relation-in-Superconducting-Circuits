import numpy as np
from matplotlib import pyplot as plt
import os

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Potentials')



def linear_potential(x, i_0):
    return i_0*x

def tilted_washboard_potential(x, i_0):
    return -i_0*x - np.cos(x)

def feynman_ratchet_potential(x, V0=1.0, lambda_=0.5, F=0.1):
    return -V0 * (np.sin(x) + (lambda_ / 4) * np.sin(2 * x)) - F * x

def sawtooth_ratchet_potential(x, V0=1.0, a=0.3, L=2*np.pi, F=0.1):

    x_mod = np.mod(x, L)  # Ensure periodicity
    V = np.zeros_like(x_mod)
    
    # Left slope (rising edge)
    mask_left = (x_mod < a * L)
    V[mask_left] = (V0 / (a * L)) * x_mod[mask_left]
    
    # Right slope (falling edge)
    mask_right = ~mask_left
    V[mask_right] = (V0 / (L - a * L)) * (L - x_mod[mask_right])
    
    # Add tilt
    V -= F * x
    
    return V

# Example usage:
x = np.linspace(0, 4 * 2*np.pi, 1000)



#############################  SAWTOOTH RATCHET  #################################
plt.figure(figsize=(10, 5))
plt.plot(x, sawtooth_ratchet_potential(x, V0=2*np.pi, a=0.5, L=2*np.pi, F=1), 'r-', linewidth=2, label="forward process, i_0 = 1")
plt.plot(x[-1]-x, sawtooth_ratchet_potential(x[-1], V0=2*np.pi, a=0.5, L=2*np.pi, F=1) -sawtooth_ratchet_potential(x, V0=2*np.pi, a=0.5, L=2*np.pi, F=1), 'b-', linewidth=2, label="backward process, i_0 = 1")

plt.plot(x, sawtooth_ratchet_potential(x, V0=2*np.pi, a=0.5, L=2*np.pi, F=0.5), 'r--', linewidth=2, label="forward process, i_0 = 0.5")
plt.plot(x[-1]-x, sawtooth_ratchet_potential(x[-1], V0=2*np.pi, a=0.5, L=2*np.pi, F=0.5) -sawtooth_ratchet_potential(x, V0=2*np.pi, a=0.5, L=2*np.pi, F=0.5), 'b--', linewidth=2, label="backward process, i_0 = 0.5")


plt.xlabel("Position (x)", fontsize=12)
plt.ylabel("Potential V(x)", fontsize=12)
plt.title("Feynman Sawtooth Ratchet Potential", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.savefig("Ratchet_Potential.pdf")
plt.show()




#############################  TILTED WASHBOARD  #################################
plt.figure(figsize=(10, 5))
plt.plot(x, tilted_washboard_potential(x, i_0 = 0.5), 'r-', linewidth=2, label="forward process, i_0 = 0.5")
plt.plot(x[-1]-x, tilted_washboard_potential(x[-1], i_0 = 0.5) -tilted_washboard_potential(x, i_0 = 0.5), 'b-', linewidth=2, label="backward process, i_0 = 0.5")

plt.plot(x, tilted_washboard_potential(x, i_0 = 1), 'r--', linewidth=2, label="forward process, i_0 = 1")
plt.plot(x[-1]-x, tilted_washboard_potential(x[-1], i_0 = 1) -tilted_washboard_potential(x, i_0 = 1), 'b--', linewidth=2, label="backward process, i_0 = 1")

plt.xlabel("Position (x)", fontsize=12)
plt.ylabel("Potential V(x)", fontsize=12)
plt.title("Tilted Washboard Potential", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.savefig("Tilted_Washboard_Potential.pdf")
plt.show()


#############################  AFFINE  #################################
plt.figure(figsize=(10, 5))
plt.plot(x, linear_potential(x, i_0 = 0.5), 'r-', linewidth=2, label="forward process, i_0 = 0.5")
plt.plot(x[-1]-x, linear_potential(x[-1], i_0 = 0.5) -linear_potential(x, i_0 = 0.5), 'b-', linewidth=2, label="backward process, i_0 = 0.5")

plt.plot(x, linear_potential(x, i_0 = 1), 'r--', linewidth=2, label="forward process, i_0 = 1")
plt.plot(x[-1]-x, linear_potential(x[-1], i_0 = 1) -linear_potential(x, i_0 = 1), 'b--', linewidth=2, label="backward process, i_0 = 1")

plt.xlabel("Position (x)", fontsize=12)
plt.ylabel("Potential V(x)", fontsize=12)
plt.title("Linear Potential", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.savefig("Linear_Potential.pdf")
plt.show()