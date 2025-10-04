import numpy as np
from matplotlib import pyplot as plt
import os

os.chdir('C:/Users/arist/Desktop/SS25/BA/Code/Fokker_Planck/FV_Solver')

def solve_fokker_planck_FV_Explicit_Euler(a_func, D: float, x: np.array, p0: np.array, T: int, dT: int, dt=None):
    # Spatial discretization
    dx = x[1] - x[0]
    N = len(p0)
    
    # Face-centered x values (i±1/2)
    x_faces = (x[:-1] + x[1:]) / 2  # length N-1
    
    # Evaluate a(x) at faces
    a_faces = a_func(x_faces)
    
    # Initialize solution array
    p = np.zeros((int(T/dT), N))
    t = np.zeros(int(T/dT))
    p[0] = p0.copy()  # Important to copy the initial condition
    p_current = p[0]
    
    # Calculate stable time step if not provided
    if dt is None:
        max_a = np.max(np.abs(a_faces))
        dt = min(0.5*dx/max_a, 0.5*dx**2/D)/2

    print('Time step size: ', dt)
    
    # Time stepping loop
    for n in range(T):
        
        # Calculate fluxes j_{i+1/2}
        p_avg = (p_current[:-1] + p_current[1:]) / 2  # length N-1
        p_upwind = np.where(a_faces > 0, p_current[:-1], p_current[1:])  #avg with upwind

        diffusion = (p_current[1:] - p_current[:-1]) * D / dx  # length N-1
        j_plus = a_faces * p_upwind - diffusion  # length N-1
        
        # Calculate fluxes j_{i-1/2}
        j_minus = np.zeros_like(j_plus)  # length N-1
        j_minus[1:] = j_plus[:-1]  # Shift right fluxes to get left fluxes
        
        # Boundary conditions (no flux)
        j_minus[0] = 0    # Left boundary
        j_plus[-1] = 0     # Right boundary
        
        # Update probability (interior points)
        p_current[1:-1] = p_current[1:-1] - (dt/dx) * (j_plus[1:] - j_minus[1:])
        
        # Boundary points (no flux)
        p_current[0] = p_current[0] - (dt/dx) * (j_plus[0] - 0)
        p_current[-1] = p_current[-1] - (dt/dx) * (0 - j_minus[-1])
        
        # Ensure non-negativity and normalization
        p_current = np.maximum(p_current, 0)
        p_current /= np.sum(p_current) * dx if np.sum(p_current) > 0 else 1

        if n%dT == 0:
            print('iteration ', n)
            p[int(n/dT)] = p_current.copy()
            t[int(n/dT)] = (n+1)*dt
    
    return p, t

# Parameters
i_0 = 1
D = 1

def a(x):
    return i_0 - np.sin(x)

xmin = -10*2*np.pi
xmax = 16000*2*np.pi
N = int((xmax-xmin)/(2*np.pi))  # ~100 points per period
N *= 5
T = 750000  # Number of time steps
dT = 5000


# Initial condition (delta function at center)
x = np.linspace(xmin, xmax, N)
p0 = np.zeros(N)
p0[int(N*(abs(xmin)/(xmax-xmin))/5)] = 1/(x[1]-x[0]) # Area = 1



# Solve
p, t = solve_fokker_planck_FV_Explicit_Euler(a, D, x, p0, T, dT)

print(p[0][N//2])
print(p[1][N//2])

# Plot
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(x, p[-100], label=rf'$\tau$={t[-100]:.2g}')
ax.plot(x, p[-50], label=rf'$\tau$={t[-50]:.2g}')
ax.plot(x, p[-1], label=rf'$\tau$={t[-1]:.2g}')
#ax.set_xlim(xmin/10, xmax/10)  # Zoom in to see details
ax.legend()
ax.grid(True)
ax.set_xlabel(r'$\varphi$')
ax.set_ylabel(r'$p(\varphi, \tau)$')
plt.tight_layout()
plt.savefig('PDF_Fokker_Planck_i_0_1_D_1.pdf')
plt.show()