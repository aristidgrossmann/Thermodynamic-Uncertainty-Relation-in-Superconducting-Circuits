import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

L = 2*np.pi
i_0 = 100
D = 0.01

x_vals = np.linspace(0, L, 100)

def P(x):
    return np.cos(x)
def U(x):
    return -i_0*x + P(x)

def lhs_integrand(y):
    return np.exp(U(y)/D)


def exponential_minus(x):
    val, _ = quad(lhs_integrand, x, x+L)
    return val

def rhs_integrand(y):
    return np.exp(i_0*y/D)


lhs = np.zeros(len(x_vals))
rhs = np.zeros(len(x_vals))

for i in range(len(x_vals)):
    lhs[i] = exponential_minus(x_vals[i])

    rhs_integral, _ = quad(rhs_integrand, x_vals[i]-L, x_vals[i])
    rhs[i] = np.exp(-P(x_vals[i]))*rhs_integral

    rhs[i] = L*np.exp(-i_0*L/(2*D))*np.exp(-i_0*x_vals[i]/D)

print(lhs)
print(rhs)
fig, ax = plt.subplots()
ax.plot(x_vals, lhs, label = 'integral', color = 'red')
ax.plot(x_vals, rhs, label = 'possible lower bound', color = 'black')
ax.legend()


plt.show()

