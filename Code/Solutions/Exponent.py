import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import i0


def f(x, y, I0, D):
    return ((np.cos(x) - np.cos(x + y) - I0 * y)/D)

def f_exp(x, y, I0, D):
    return np.exp((np.cos(x) - np.cos(x + y) - I0 * y)/D)

def integrand(x, i_0, D):
        return np.exp(-i_0/D*x)*i0(2*np.sin(x/2)/D)
L = 2*np.pi

# Create a grid of (x, y) values
x = np.linspace(0, L, 100)
y = np.linspace(0, L, 100)
X, Y = np.meshgrid(x, y)

i_0 = 2  # Example parameter value
D = 0.005
Z = f_exp(X, Y, i_0, D)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)')
ax.set_title(f'$f(x, y) = \cos(x) - \cos(x+y) - {i_0}x$')

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()


fig, ax = plt.subplots()
ax.plot(x/np.pi, integrand(x, i_0, D), color = 'red')
ax.plot(x/np.pi, D*(1-np.exp(-i_0*2*np.pi/D))*np.ones_like(x), color = 'blue')
#ax.set_yscale('log')
ax.set_xlabel(r'x/\pi')
ax.grid(True)
plt.show()


