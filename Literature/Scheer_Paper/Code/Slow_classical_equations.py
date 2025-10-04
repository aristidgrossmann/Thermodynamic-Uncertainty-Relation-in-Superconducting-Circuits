#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:08:50 2024

@author: david
"""

"""
This code integrates the slow classical equations of motion for the effective degrees of freedom of 
the clock circuit. It samples the relevant Langevin equation for many realizations of the nyquist Johnson 
noise to obtain the uncertainty product as well as the variance of the light phase and the 
mean amplitude of the oscillation as well as the resulting DC voltage for the counter.
"""
import numpy as np
from numba import njit
from tan_solver import tan_solution
import matplotlib.pyplot as plt

#Set up random number generator
rng = np.random.default_rng(0)

#%%
# =============================================================================
# Parameters of the Simulation
# =============================================================================
dt=0.5 #Size of the timesteps

#Parameters from classical equation
gamma = 0.1 # Damping of the oscillator
Gamma = 0.001 # Damping of the counter
kappa = 0.01  #coupling strength (corresponds to Hamiltonian timescale divided by large factor)

epsilon=0.0005 #noise strength

#%%
# =============================================================================
# Auxiliary parameters
# =============================================================================
w = 1 # Resonance frequency
f = Gamma*w # External force (current) matching the frequency for f=Gamma*w
dw=1-f/Gamma #Frequency mismatch
y=gamma
a=epsilon
c=a

# Get the the steady state solution for synchronization without noise
alpha=2*(dw)/y
beta=kappa**2/y**2/Gamma

t0=(tan_solution(alpha,beta)[0])[0]

#%%
# =============================================================================
# Functions
# =============================================================================
# A single timestep with the forward euler method
@njit
def timestep(theta,phi,A,df,y=y,Gamma=Gamma,dw=dw,kappa=kappa):
    theta_new=theta+dt*(-dw+kappa*A*np.cos(phi-theta)/(2*Gamma))+df
    phi_new=phi+dt*kappa*np.sin(phi-theta)/(2*A)
    A_new=A+dt*(-y*A/2-kappa*np.cos(phi-theta)/2)
    return theta_new,phi_new,A_new

# Integration with the forward Euler method for a given noise configuration
@njit 
def integrate(times,theta_0,phi_0,A0,f_noise,y=y,Gamma=Gamma,dw=dw,kappa=kappa):
    theta_vals=np.zeros(len(times),dtype=np.float64)
    phivals=np.zeros(len(times),dtype=np.float64)
    Avals=np.zeros(len(times),dtype=np.float64)
    theta,phi,A=theta_0,phi_0,A0
    for j,df in enumerate(f_noise):
        theta,phi,A=timestep(theta,phi,A,df,y=y,Gamma=Gamma,dw=dw,kappa=kappa)

        theta_vals[j]=theta
        phivals[j]=phi
        Avals[j]=A
    return theta_vals-theta_0,phivals-phi_0,Avals #subtract theta_0 to calculate the average slope

#Calculate the steady state without noise to use as an initial condition
@njit
def initials(t0,times,f_noise,y=y,Gamma=Gamma,dw=dw,kappa=kappa):
    theta0=0
    phi0=np.arctan(t0) #set phi to 0 since only the phase difference is fixed
    A0=-kappa/(y*np.sqrt(1+t0**2))
    theta,phi, A=integrate(times,theta0,phi0,A0,f_noise,y=y,Gamma=Gamma,dw=dw,kappa=kappa)
    return theta0,phi0,A[-1]


#Obtain mean and variance over time traces
@njit
def trace_stats(tr_lists,times):
    l=len(tr_lists)
    mean=np.sum(tr_lists,axis=0)/l
    var=np.sum(tr_lists**2,axis=0)/l-mean**2
    return mean,var

# Generate time-traces for theta,phi and A for different noise configurations
def Sample(y=y,Gamma=Gamma,dw=dw,kappa=kappa,a=a,tmax=6000,nsamples=4000,lbar=True):
    #Set up time
    times=np.arange(0,tmax,dt)
    #Generate independent noise for all time traces
    noise_traces=rng.standard_normal(size=(nsamples,len(times)))*np.sqrt(dt)*np.sqrt(a)
    
    #Set up arrays for the different time traces
    theta_traces=np.empty((nsamples,len(times)))
    phi_traces=np.empty((nsamples,len(times)))
    A_traces=np.empty((nsamples,len(times)))
    
    #Calculate initial conditions from the expected steady state without noise
    alpha=2*(dw)/y
    beta=kappa**2/y**2/Gamma
    t0=(tan_solution(alpha,beta)[0])[0] #Select first stable solution since it usually corresponds to the case of frequency locking
    
    #Simulate with good initial conditions to reach steady state
    init_noise=rng.standard_normal(len(times))*np.sqrt(dt)*np.sqrt(a)
    theta0,phi0,A0=initials(t0,times,init_noise,y,Gamma,dw,kappa)
    if lbar:
        for i,noise in enumerate(noise_traces):
            theta_traces[i],phi_traces[i],A_traces[i]=integrate(times,theta0,phi0,A0,noise,y=y,Gamma=Gamma,dw=dw,kappa=kappa)
    else:
        for i,noise in noise_traces:
            theta_traces[i],phi_traces[i],A_traces[i]=integrate(times,theta0,phi0,A0,noise,y=y,Gamma=Gamma,dw=dw,kappa=kappa)
    return theta_traces,phi_traces,A_traces,times

#Obtain the values of the uncertainty product and the variance of the oscillator phase along th IV-curve
def f_sweep(fmin=0.6,fmax=1.6,nsteps=500,y=y,Gamma=Gamma,kappa=kappa,a=a,tmax=5000,nsamples=80000,lbar=False):
    fvals=np.linspace(fmin,fmax,nsteps)
    dw_vals=1-fvals/Gamma
    
    TUR_vals=np.zeros(len(fvals))
    U_vals=np.zeros(len(fvals))
    phivar_vals=np.zeros(len(fvals))
    var_vals=np.zeros(len(fvals))
    
    for i,dw in enumerate(dw_vals):
        theta_traces,phi_traces,A_traces,times=Sample(y=y,Gamma=Gamma,dw=dw,kappa=kappa,a=a,tmax=tmax,nsamples=nsamples,lbar=lbar)
        theta_mean=np.mean(theta_traces[::,-1])
        theta_var=np.var(theta_traces[::,-1])
        phi_var=np.var(phi_traces[::,-1])

        var_vals[i]=theta_var/tmax
        phivar_vals[i]=phi_var/tmax
        TUR_vals[i]=2*fvals[i]/Gamma/a*theta_var/(times[-1]+theta_mean)
        U_vals[i]=1+theta_mean/tmax
    
    return TUR_vals,U_vals,fvals,phivar_vals,var_vals


#%%
# =============================================================================
# Simulate results for a sweep of the external current 
# =============================================================================

TUR,U_vals,f,varphi,var=f_sweep(fmin=0.01*Gamma,fmax=1.7*Gamma)

#%%
#Save simulation data
# np.save('2ndTUR_withQ_y_0.1_k_0.01_G_0.001_T_0.0005_steps_500_tmax_5000_nsamp_80000',np.array([f,U_vals,TUR,varphi]))

#Plotting
fig, axes = plt.subplots(nrows=3,sharex=True)
axes[0].plot(f,U_vals)
axes[0].plot(f,f/Gamma,linestyle='--')

axes[1].plot(f,TUR)
axes[1].plot(f,2*np.ones(len(f)),linestyle='--')

axes[1].plot(f,4*y/beta*np.ones(len(f)),color='black',linestyle='dotted')

ax2=axes[1].twinx()
ax2.plot(f,np.gradient(U_vals,f)**2*Gamma**2,color='#2ca02c',linestyle='-.')

axes[1].set_ylim(0,2.2)
ax2.set_ylim(0,1.1)
ax2.tick_params(axis='y', colors='#2ca02c')

axes[2].plot(f,U_vals/varphi)
axes[2].plot(f,U_vals/(c*np.ones(len(f))),linestyle='--')
axes[2].plot(f,1/(c*np.ones(len(f))/(1+beta/2)**2),linestyle='dotted',color='black')

axes[2].set_yscale('log')

plt.show()




