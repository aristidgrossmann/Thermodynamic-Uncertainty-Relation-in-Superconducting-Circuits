#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:41:05 2024

@author: david
"""


"""
This code solves the coupled Lindblad equation and Langevin equation for a given
set of System parameters (including temperature) to obtain an IV-curve of the counter.
It compares the result to the analytical result in the classical limit to show,
when it becomes valid.
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from numpy.random import default_rng
from scipy.special import hyp1f1
from tan_solver import calculate_curves 

#Set up random number generator
rng = default_rng(0)

#%%
# =============================================================================
# Parameters of the simulation
# =============================================================================
#
M = 200 # Size of Hilbertspace (adapt to your expected photon number on resonance!)

dt = 0.01 #Size of the timesteps

#Parameters from classical equation
gamma = 0.1 # Damping of the oscillator
Gamma = 0.001 # Damping of the counter
kappa = 0.01  #coupling strength (corresponds to Hamiltonian timescale divided by large factor)

r_vals=[0.01,0.02,0.05]  #Square root of different values of the Light matter coupling (Note that this the square root of the r used in the paper)

#Thermal occupation of the Oscillator (Set to zero in our simulations)
n0 = 0 


#%%
# =============================================================================
# Auxiliary parameters
# =============================================================================
w = 1 # Resonance frequency (set to one by choosing the timescale)
r = 0.01
f = Gamma*w # External force (current) matching the frequency for f=Gamma*w
i = 1j  # imaginary unit
# System parameters
G = Gamma/(2*w*r**2) #ratio between damping constant of counter and quantum conductance

dw = w-f/Gamma # Frequency mismatch (tunable by external current)
h = kappa/(2*w*r**2) #Hamiltonian timescale

T = 1 / np.log(1+1/n0) if n0 != 0 else 0 #kbT/hbar w
#Noise strength for Langevin equation
c=np.sqrt(Gamma*T)/r

#%%
# =============================================================================
# Functions
# =============================================================================
#matrices
a_diag = np.sqrt(np.linspace(1,M-1,M-1), dtype = complex)
a = sp.csc_matrix(sp.diags(a_diag, 1))
ad = sp.csc_matrix(sp.diags(a_diag, -1))
one = sp.csc_matrix(sp.identity(M))


#Set up number operator explicitly to avoid rounding errors
n_diag = np.array(np.linspace(0,M-1,M),dtype =complex)
n= sp.csc_matrix(sp.diags(n_diag))

#Define operator valued Bessel function (not normal ordered yet)
def Bessel(r=r):
    diag=r*hyp1f1(-np.linspace(0,M-1,M),2,r**2) #Normal ordered Bessel function
    return sp.csc_matrix(sp.diags(diag))

print(Bessel(r))

# Calculate the density matrix in the next timestep using the Lindblad equation
def Lindblad(rho,theta, dt,r=r,gamma=gamma,h=h):
    #Define Hamiltonian
    H = i*h*(Bessel(r)@a*np.exp(i*np.real(theta))-ad@Bessel(r)*np.exp(-i*np.real(theta)))*np.exp(-r**2/2)/2
    #Lindblad timestep 
    rho_new =  dt*(-i*(H@rho - rho@H) + gamma*(n0+1)*(a@rho@ad - (ad@a@rho + rho@ad@a)/2)+gamma*n0*(ad@rho@a - (a@ad@rho + rho@a@ad)/2)) +rho
   
    return rho_new

#Calculate theta in the next timestep using the Langevin equation
def Langevin(rho,theta, dt,c=c,r=r,G=G,dw=dw):
    O=Bessel(r)@a*np.exp(i*theta)+ad@Bessel(r)*np.exp(-i*theta)
    theta_new=theta+dt*(-G*dw+h/2*np.exp(-r**2/2)*np.real((O@rho).trace()))/G+np.sqrt(dt)*c*rng.standard_normal(1)[0]/G

    return theta_new

def initialize(r=r,gamma=gamma,G=G,dw=dw,h=h,c=c):
    tmax = 50
    times = np.arange(0,tmax,dt)
    #initial state=vacuum state
    init_diag = np.zeros(M)
    init_diag[0]=1
    init = sp.csc_matrix(sp.diags(init_diag, 0))
    #set initial conditions
    rho_old=init
    theta_old=0

    #time evolution
    for j in range(len(times)):
        rho_new = Lindblad( rho_old, theta_old, dt,r=r,gamma=gamma,h=h)
        theta_new=Langevin(rho_old,theta_old,dt,c=c,r=r,G=G,dw=dw)
        rho_old = rho_new
        theta_old=theta_new
        
    return rho_new,theta_new


def evolve_system(rho,theta,tmax=40,r=r,gamma=gamma,G=G,dw=dw,h=h,c=c):
    tmax = tmax
    times = np.arange(0,tmax,dt)
    rho_old=rho
    theta_old=theta
    nvals = np.zeros(len(times), dtype= complex)
    
    for j in range(len(times)):
        rho_new = Lindblad( rho_old, theta_old, dt,r=r,gamma=gamma,h=h)
        theta_new=Langevin(rho_old,theta_old,dt,c=c,r=r,G=G,dw=dw)
        nvals[j]=(n@rho).trace()
        rho_old = rho_new
        theta_old=theta_new
        
    return rho_new,theta_new,np.mean(nvals)


#Simulate the resulting IV curve with potential hysteresis for different values of r
def IV_curve(fmin=0,fmax=6*w*Gamma,fsteps=120,r_vals=[0.01,0.02,0.05],tmax=40):
    #Set up f_vals to allow for hysteresis
    vals=np.linspace(fmin,fmax,fsteps)
    f_vals=np.concatenate((vals,vals[-2::-1]),axis=None)
    
    nmax=np.empty(len(r_vals))
    voltage_vals=np.empty((len(r_vals),len(f_vals)))
    photon_numbers=np.empty((len(r_vals),len(f_vals)))
    
    for i,r in enumerate(r_vals):
        G=Gamma/(2*w*r**2) #ratio between damping constant of counter and quantum conductance
        h=kappa/(2*w*r**2) #Hamiltonian timescale
        c=np.sqrt(Gamma*T)/r #Noise strength
        #Simulate an initial state
        nmean=n0+kappa**2/(4*r**2*w**2*gamma**2) #expected photon number in the classical regime

        print('maximum mean photon number: ',str(nmean))
        print('r^2n: ',str(r**2*nmean))
        nmax[i]=nmean
        rho_0,theta_0=initialize(r ,gamma , G, w, h ,c)
        
        for j,f in enumerate(f_vals):
            dw=w-f/Gamma#frequency mismatch (tunable by external current)
            rho,theta,n_val=evolve_system(rho_0,theta_0,tmax=tmax,r=r,gamma=gamma,G=G,dw=dw,h=h,c=c)
            voltage_vals[i,j]=(theta-theta_0)/tmax
            photon_numbers[i,j]=np.real(n_val)
            rho_0=rho
            theta_0=theta

    return r_vals,f_vals,voltage_vals,photon_numbers,nmax
#%%
#Run IV curve simulation
r_vals,f_vals,voltage_vals,photon_numbers,nmax=IV_curve(fmin=0*w*Gamma,fmax=6*w*Gamma,fsteps=120,r_vals=r_vals,tmax=40)
#%%
#Save file
# np.save('IV_curve_y_0.1_k_0.01_G_0.001_T_0_r_0.0(125)_steps_120_tmax_40',np.array([np.array([f_vals,f_vals,f_vals]),voltage_vals,photon_numbers]))

#Plotting
vals=np.linspace(0,6*w*Gamma,120)
f_vals=np.concatenate((vals,vals[-2::-1]),axis=None)

#Obtain analytical results
fvals_2=np.linspace(0,6*w*Gamma,100000)
avals=2*(w-fvals_2/Gamma)/gamma
beta=kappa**2/gamma**2/Gamma/w
tvals,avals=calculate_curves(avals, beta)
new_fvals=[]
for a in avals:
    new_fvals.append((w-gamma*a/2)*Gamma)
    
plt.figure()
for i,r in enumerate(r_vals):
    plt.plot(f_vals,voltage_vals[i]+w)
plt.legend(r_vals) 
for i in range(3):
    plt.plot(new_fvals[i],w-gamma/2*tvals[i])
  
