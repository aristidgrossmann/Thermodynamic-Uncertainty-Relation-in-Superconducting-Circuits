import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import i0
from matplotlib import pyplot as plt


#########################   HELPER FUNCTIONS    ###################################
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
    result, _ = quad(I_tilde_plus_minus, 0, L, args= (i_0, D))
    return result/(D*(1 - np.exp(- i_0*L/D)))

#####  FPT Variance  ########
def Delta_T_2_paper(i_0, D):
    L = 2*np.pi
    
    def integrand_z(z, i_0, D):
        return I_tilde_plus_minus(z, i_0, D)**2 * I_tilde_minus_plus(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return 2/D**2 /(1 - np.exp(- i_0*L/D))**3 * result 

def Delta_T_2_me(i_0, D):
    L = 2*np.pi
    
    def integrand_alpha(alpha, i_0, D):
        result = I_tilde_plus_minus(alpha, i_0, D)**2 * np.exp(-V(alpha, i_0, D))
        return result
    
    def integrand_phi(phi, i_0, D):
        result, _ = quad(integrand_alpha, phi - L, phi, args=(i_0, D))
        return result*np.exp(V(phi, i_0, D))

    result, _= quad(integrand_phi, 0, L, args= (i_0, D))
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
        return I_tilde_plus_minus(z, i_0, D) * I_tilde_minus_plus(z, i_0, D)

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return L/(D*(1 - np.exp(- i_0*L/D)))**2 * result /T_1_paper(i_0, D)**2

def derivative_r_d(i_0, D):

    def inner_integrand(x, i_0, D):
        return np.exp(-V(x, i_0, D))*I_tilde_plus_minus(x, i_0, D)
    
    def outer_integrand(y, i_0, D):
        inner_integral, _ = quad(inner_integrand, y-L, y, args = (i_0, D))

        return np.exp(V(y, i_0, D))*I_tilde_minus_plus(y, i_0, D)*inner_integral
    
    outer_integral, _ = quad(outer_integrand, 0, L, args = (i_0, D))
    first_term = -L*2.0/D/(1-np.exp(-i_0*L/D))**2*outer_integral/T_1_paper(i_0, D)**2
    second_term = 2.0*T_1_paper(i_0, D)/L *r_d(i_0, D)**2

    return  second_term
    


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


#########################   testing   ###################################
def Unc_prod(i_0, D):

    def integrand_z(z, i_0, D):
        return 1.0/I_tilde_minus_plus(z, i_0, D)**2

    result, _= quad(integrand_z, 0, L, args= (i_0, D))
    return 2*L*i_0/((D*(1 - np.exp(- i_0*L/D)))**2 * result )







#########################   APPROX   ###################################
def compute_V_D_larger_than_2(i_0, D):
    return i_0/(1+ D/(2*(i_0**2 + D**2)))

def compute_V_i0_larger_than_D(i_0, D):
    L = 2*np.pi
    return i_0*(1-np.exp(-L*i_0/D))/(1+ 1/(2*i_0**2))




if __name__ == 'main':
    B = 10
    D = B**2 /2
    D = 0.1
    i_0_vals = np.linspace(0, 2, 11)




    V_vals_double_integral = np.zeros_like(i_0_vals)
    for i in range(len(i_0_vals)):
        V_vals_double_integral[i] = L*Probability_Current_Bessel(i_0_vals[i], D)

    V_vals_Bessel = np.zeros_like(i_0_vals)
    for i in range(len(i_0_vals)):
        V_vals_Bessel[i] = L*Probability_Current_Bessel(i_0_vals[i], D)


    ## approximations
    # V_vals_D_larger_than_2 = np.zeros_like(i_0_vals)
    # for i in range(len(i_0_vals)):
    #     V_vals_D_larger_than_2[i] = compute_V_D_larger_than_2(i_0_vals[i], D)

    # V_vals__i0_larger_than_D = np.zeros_like(i_0_vals)
    # for i in range(len(i_0_vals)):
    #     V_vals__i0_larger_than_D[i] = compute_V_i0_larger_than_D(i_0_vals[i], D)



    fig, ax = plt.subplots()

    #ax.plot(i_0_vals, analytical_Voltage_mean_no_Noise(i_0_vals), label = 'no Noise', color = 'blue')
    # ax.scatter(i_0_vals, V_vals_Bessel_larger_exponent, marker = 'x', s = 10, label = 'Noise Bessel larger exponent', color = 'red')
    ax.scatter(i_0_vals, V_vals_double_integral, marker = 'x', s = 10, label = 'Noise double integral', color = 'red')
    ax.plot(i_0_vals, V_vals_Bessel, label = 'Noise Bessel', color = 'black')

    ## approximations
    # ax.scatter(i_0_vals, V_vals_D_larger_than_2, marker = 'x', s = 10, label = 'Noise D > 2 approx', color = 'blue')
    # ax.scatter(i_0_vals, V_vals__i0_larger_than_D, marker = 'x', s = 10, label = 'Noise i_0 >> D approx', color = 'purple')



    ax.grid(True)
    ax.legend()
    plt.show()

