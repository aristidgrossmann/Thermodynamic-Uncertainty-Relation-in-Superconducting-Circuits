import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import i0
from matplotlib import pyplot as plt


#########################   HELPER FUNCTIONS    ###################################
L = 1

def V(phi, i_0, D):
    return (-i_0*phi - np.cos(phi))/D

def I_plus_minus(x, i_0, D):

    def integrand(alpha, i_0, D):
        return np.exp(- V(alpha, i_0, D))
    
    integral, _ = quad(integrand, x - L, x, args= (i_0, D))
    return integral*np.exp( V(x, i_0, D))

def I_minus_plus(x, i_0, D):

    def integrand(alpha, i_0, D):
        return np.exp( V(alpha, i_0, D))
    
    integral, _ = quad(integrand, x , x + L, args= (i_0, D))
    return integral*np.exp( -V(x, i_0, D))

def I_tilde_plus_minus(x, i_0, D):
    return 1/(D*(1-np.exp(-i_0*L/D)))*I_plus_minus(x, i_0, D)

def I_tilde_minus_plus(x, i_0, D):
    return 1/(D*(1-np.exp(-i_0*L/D)))*I_minus_plus(x, i_0, D)



#########################   Ratchet Analytical solution helper function    ###################################
V_0 = 2
a = 0.3
b = L-a

def V_Ratchet(x, i_0, D):
    ratchet_potential = 0
    z = x % L
    if z <= a:
        ratchet_potential = V_0*z/a
    else:
        ratchet_potential = V_0*(1 - (z-a)/b)
    
    return (-i_0*x + ratchet_potential)/D

### V_tilde
def V_tilde_a(x, i_0):
    return x*(i_0 - V_0/a)

def V_tilde_b(x, i_0):
    return x*(i_0 + V_0/b)

### E_a,b_plus_minus
def E_a_plus(x, i_0, D):
    return np.exp(V_tilde_a(x, i_0)/D)

def E_a_minus(x, i_0, D):
    return np.exp(-V_tilde_a(x, i_0)/D)

def E_b_plus(x, i_0, D):
    return np.exp(V_tilde_b(x, i_0)/D)

def E_b_minus(x, i_0, D):
    return np.exp(-V_tilde_b(x, i_0)/D)

def I_plus_minus_Ratchet(x, i_0, D):
    A_a_minus = (1-np.exp(-V_0*L/(D*a))*E_a_minus(b, i_0, D))/V_tilde_a(1, i_0) + (E_b_minus(b, i_0, D)-1)/V_tilde_b(1, i_0)
    B_a_minus = (np.exp(-V_0*L/(D*a))*E_a_minus(L, i_0, D)-1)/V_tilde_a(1, i_0)

    A_b_plus = np.exp(V_0/D*(1+a/b))*(1- E_a_plus(a, i_0, D))/V_tilde_a(1, i_0) + (E_b_plus(a, i_0, D)-np.exp(V_0*L/(D*b)))/V_tilde_b(1, i_0) 
    B_b_plus = (np.exp(V_0*L/(D*b))*E_b_plus(-L, i_0, D)-1)/V_tilde_b(1, i_0)


    if x <= a:
        result =  - D*(E_a_minus(x, i_0, D)*A_a_minus + B_a_minus)
    else: 
        result =  -D*(E_b_minus(x, i_0, D)*A_b_plus + B_b_plus)
    return result

def I_minus_plus_Ratchet(x, i_0, D):
    A_a_plus = (1-np.exp(V_0*L/(D*a))*E_a_plus(b, i_0, D))/V_tilde_a(1, i_0) + (E_b_plus(b, i_0, D)-1)/V_tilde_b(1, i_0)
    B_a_plus = (np.exp(V_0*L/(D*a))*E_a_plus(L, i_0, D)-1)/V_tilde_a(1, i_0)

    A_b_minus = np.exp(-V_0/D*(1+a/b))*(1- E_a_minus(a, i_0, D))/V_tilde_a(1, i_0) + (E_b_minus(a, i_0, D)-np.exp(-V_0*L/(D*b)))/V_tilde_b(1, i_0) 
    B_b_minus = (np.exp(-V_0*L/(D*b))*E_b_minus(-L, i_0, D)-1)/V_tilde_b(1, i_0)

    
    if x <= a:
        result =   D*(E_a_plus(x, i_0, D)*A_a_plus + B_a_plus)
    else: 
        result =  D*(E_b_plus(x, i_0, D)*A_b_minus + B_b_minus)
    
    return result*np.exp(- i_0*L/D)
    
def T_1_Ratchet(i_0, D):
    result, _ = quad(I_minus_plus_Ratchet, 0, L, args = (i_0, D))
    return result/(D*(1 - np.exp(- i_0*L/D)))

def Delta_T_2_Ratchet(i_0, D):
    
    def integrand_z(z, i_0, D):
        return I_plus_minus_Ratchet(z, i_0, D)**2 * I_minus_plus_Ratchet(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return 2/D**2 /(1 - np.exp(- i_0*L/D))**3 * result 

def Uncertainty_Product_Ratchet(i_0, D):
    return L*Delta_T_2_Ratchet(i_0, D) / T_1_Ratchet(i_0, D)**2 * i_0/D




#########################   Probability Current   ###################################
def Probability_Current_double_integral(i_0, D):

    def integrand(x, y):
        return np.exp((np.cos(x)-np.cos(x+y) - i_0*y)/D)
    
    integral, _ = dblquad(integrand,0, L, 0, L)
   
    return D*(1-np.exp(-L*i_0/D))/integral

def analytical_Voltage_MFPT_double_integral(i_0, D):

    def integrand(x, y):
        return np.exp((-np.cos(x)+np.cos(x-y) - i_0*y)/D)
    
    integral, _ = dblquad(integrand,0, L, 0, L)
    denominator = L*D*(1-np.exp(-L*i_0/D))
    V = denominator/integral    
    return V


def Probability_Current_Bessel(i_0, D):

    def integrand(x):
        return np.exp(-i_0/D*x)*i0(2*np.sin(x/2)/D)
    
    integral, _ = quad(integrand,0, L)
    denominator = D*(1-np.exp(-L*i_0/D))
    V = denominator/(integral*L)
    return V




#########################   FPT   ###################################
#####  FPT Mean  ########
def T_1_paper(i_0, D):
    result, _ = quad(I_plus_minus, 0, L, args= (i_0, D))
    return result/(D*(1 - np.exp(- i_0*L/D)))

#####  FPT Variance  ########
def Delta_T_2_paper(i_0, D):
    
    def integrand_z(z, i_0, D):
        return I_plus_minus(z, i_0, D)**2 * I_minus_plus(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return 2/D**2 /(1 - np.exp(- i_0*L/D))**3 * result 


#########################   Phase Cumulant Velocities   ###################################
#####  Phase Mean Velocity ########
def Mean_Velocity(i_0, D):
    return L /T_1_paper(i_0, D)

#####  Phase Variance Velocity ########
def Variance_Velocity(i_0, D):
    return  L**2*Delta_T_2_paper(i_0, D) / T_1_paper(i_0, D)**3

#########################   Uncertainty Product   ###################################
def Uncertainty_product(i_0, D):
    return L*Delta_T_2_paper(i_0, D) / T_1_paper(i_0, D)**2 * i_0/D


#########################   differential resistance   ###################################
def r_d(i_0, D):

    def integrand_z(z, i_0, D):
        return I_plus_minus(z, i_0, D) * I_minus_plus(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return L/(D*(1 - np.exp(- i_0*L/D)))**2 * result /T_1_paper(i_0, D)**2
    

#########################   analytical voltage, no noise    ###################################
def analytical_Voltage_mean_no_Noise(i):
    vals = np.zeros_like(i)
    for j in range(len(i)): 
        if i[j] <= 1:
            continue
        else: 
            vals[j] = np.sqrt(i[j]**2 -1)
    #return np.where(i <= 1, 0, np.sqrt(i**2 -1))
    return vals


