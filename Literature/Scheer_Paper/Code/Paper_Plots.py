#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:01:02 2024
@author: david
"""

"""
This code generates the plots unsed in our paper up to some cosmetics in inkscsape
"""
import numpy as np
import matplotlib.pyplot as plt
from tan_solver import calculate_curves

#Simulated parameters for the paper
Gamma=0.001  #Damping of the counter
gamma=0.1    #Damping of the oscillator
kappa=0.01   #Coupling constant
w=1          #Resonance frequency

#%%
# =============================================================================
# Plot for Figure 2 of the Paper (IV-curves)
# =============================================================================
r_vals=[0.01,0.02,0.05]
f_vals,voltage_vals,photon_numbers = np.load('IV_curve_y_0.1_k_0.01_G_0.001_T_0_r_0.0(125)_steps_120_tmax_40.npy')

#Plotting
f_vals=f_vals[0]/Gamma #extract one of the f-arrays since they are all the same

#Obtain analytical results
fvals_2 = np.linspace(0,6*w*Gamma,100000)
avals = 2*(w-fvals_2/Gamma)/gamma
beta = kappa**2/gamma**2/Gamma/w
tvals,avals = calculate_curves(avals, beta)
new_fvals = []

for a in avals:
    new_fvals.append((w-gamma*a/2))
    
plt.figure()
linestyles=['-','--','-.']
for i,r in enumerate(r_vals):
    plt.plot(f_vals,voltage_vals[i]+w,linestyle=linestyles[i])
plt.legend([0.0001,0.0004,0.0025])

#Plot analytical results
colors=['black','black','grey'] #Set up color scheme for analytical results
for i in range(3):
    plt.plot(new_fvals[i],w-gamma/2*tvals[i],color=colors[i],linestyle='--')
   
plt.xlabel(r'$f/(\Gamma\omega_0)$') 
plt.ylabel(r'$2eV_0/(\hbar\omega_0)$')   
plt.xlim(0,4)
plt.ylim(0,5)

plt.show()

#%%
# =============================================================================
# Plot for Figure 3 of the Paper (TUR-violation)
# =============================================================================
#Reload file with simulation data
f,U_vals,TUR,varphi=np.load('TUR_withQ_y_0.1_k_0.01_G_0.001_T_0.0005_steps_500_tmax_5000_nsamp_80000.npy')

epsilon=0.0005  #Noise strength determined by the Temperature

beta=kappa**2/(2*Gamma*gamma**2) #Synchronization strength

f=f/Gamma #rescale x-Axis

fig, axes = plt.subplots(nrows=3,sharex=True)
#Plot IV-curve
axes[0].plot(f,U_vals)  #IV-curve of the clock circuit
axes[0].plot(f,f,linestyle='--')  #IV-curve of an Ohmic resistance
axes[0].set_ylabel(r'$2eV_0/(\hbar\omega_0)$')

#Plot Uncertainty product and squared differential resistance
axes[1].plot(f,TUR)
axes[1].plot(f,2*np.ones(len(f)),linestyle='--')
axes[1].plot(f,4*gamma/beta*np.ones(len(f)),color='black',linestyle='dotted')
axes[1].set_ylabel(r'$\langle\!\langle Y^2\rangle\!\rangle\sigma t/(k_B\langle Y\rangle^2)$')
axes[1].set_ylim(0,2.2)

ax2=axes[1].twinx() #Twin axis for differential resistance

ax2.plot(f,np.gradient(U_vals,f)**2,color='#2ca02c',linestyle='-.')
ax2.set_ylim(0,1.1)
ax2.tick_params(axis='y', colors='#2ca02c')
ax2.set_ylabel(r'$(GR_Q)^2$',color='#2ca02c')

#Plot resulting quality facotr
axes[2].plot(f,U_vals/varphi)

axes[2].plot(f,1/(epsilon*np.ones(len(f))),linestyle='--')
axes[2].plot(f,1/(epsilon*np.ones(len(f))/(1+beta)**2),linestyle='dotted',color='black')

axes[2].set_yscale('log')
axes[2].set_ylim(10e2,13e4) #cut off the off-resonant quality factors for better visibility of the resonance region

axes[2].set_ylabel(r'$Q$')
axes[2].set_xlabel(r'$f/(\Gamma\omega_0)$') 

plt.show()

