
"""
Created on Thu Jun 20 14:47:56 2024

@author: david
"""

# Simulation code for our paper "The superconducting clock-circuit: improving
 the coherence of Josephson radiation beyond the thermodynamic uncertainty relation"
A Python project for simulating the superconducting clock circuit and comparing
to the analytical results in the classical limit.


## Introduction
This project simulates the superconducting clock circuit 
using either a Lindblad equation for the oscillator coupled to a Langevin
equation for the counter or a set of Langevin equations in the classical limit.
 

## Features
- Simulate quantum and classical equations of motion with custom circuit parameters.
- Calculate results to noiseless analytics in the classical limit.
- Visualize simulation results.


## Installation
### Prerequisites
Essential:
- Python 3.8+ (We used version 3.11.5)
- NumPy (version 1.24.3)
- SciPy (version 1.11.1)
- Matplotlib (version 3.7.2)

Optional:
- numba (version 0.57.1)
# This speeds up the simulation considerably. If not wanted, 
# simply remove the @njit wrappers from all functions


## Usage
###Files
#------------------------------------------
    "tan_solver.py"
    This code solves the cubic equation (t^2+a)(t+1)=b given by the synchronization condition.
    In addition, it checks the stability of the obtained real solutions in order to obtain
     a full picture of the steady state.
     
       Functions:
       - "tan_soultion" solves the equation for a given set of parameters (a,b)
         and returns all real solutions classified by their stability.

       - "calculate_curves" takes a list of a-values for a given b and returns the 
         categorized real solutions as a function of a.
#------------------------------------------
    "recover_classical_equations.py"
    This code solves the coupled Lindblad equation and Langevin equation for a given
    set of System parameters (including temperature) to obtain an IV-curve of the counter.
    It compares the result to the analytical result in the classical limit to show,
    when it becomes valid. In the end, it plots the simulation results to compare them to the 
    analytics in the classical limit. The results with the default parameters correspond to 
    Figure 2 of the paper.
    
        Input parameters:
        -M (Size of the Hilbert space. Adapt to fit your expected photon number+variance!)
        -dt (Size of discrete timestep. Adapt, when dealing with faster timescales!)
        
        -gamma  (damping constant of the oscillator)
        -Gamma  (damping constant of the counter)
        -kappa  (coupling constant)
        -r_vals (list of values of the SQUARE ROOT!! of the light-matter coupling constant)
    
        -n0     (Thermal occupation of the oscillator. Is set to zero in our paper.)
    
         Functions:
         - "initialize" integrates the equations of motion for a time tmax starting from the vacuum state.
           This function is supposed to initialize the steady state for a given set of parameters

         - "evolve_system" integrates an arbitrary input state (rho,theta) for a time tmax. 
           It takes parameters that are rescaled to give constant ceffective parameters for varying
           light-matter coupling r^2.
         
         -"IV_curve" integrates the euqations of motion with a given set of r-values for varying bias current f. 
          It always takes the last result as the initial step for the next one and sweeps f up and then down to detect hysteresis.
          For each value of r, the function rescales the parameters of the equations of motin in order to maintain
           the effective input parameters at a constant value. For each r, it prints the expected photon number at resonance.
           This indicates, how large the Hilbert space needs to be for the simulation.
#------------------------------------------
    "Slow_classical_equations.py"
    This code integrates the slow classical equations of motion for the effective degrees of freedom of 
    the clock circuit. It samples the relevant Langevin equation for many realizations of the nyquist Johnson 
    noise to obtain the uncertainty product as well as the variance of the light phase and the 
    mean amplitude of the oscillation as well as the resulting DC voltage for the counter. The results 
    with the default parameters correspond to Figure 3 of the paper.
    
    Input parameters:
    -dt (Size of discrete timestep. Adapt, when dealing with faster timescales!)
    
    -gamma   (damping constant of the oscillator)
    -Gamma   (damping constant of the counter)
    -kappa   (coupling constant)
    -epsilon (strength of Nyquist-Johnson noise)

    Functions:
    - "integrate" integrates the equations of motion to obtain a time trace
      for a given noise configuration. For theta and phi,
      it outputs the change with respect to the initial value
      
    - "initials" calculates initial conditions that correspond to the expected steady state
    
    - "trace_stats" calculates the mean and variance in each timestep for a given list of 
      time traces of a variable.
      
    - "Sample" simulates a list of time traces for a desired number nsample of noise realizations.
    
    - "f_sweep" simulates a full IV-curve in a given range of the bias current f. It also calculates
      the uncertainty product as well as the variances of phi and theta at each point.
#------------------------------------------
    "Paper_Plots.py"    
    This code generates the plots unsed in our paper up to some cosmetics in inkscsape.
    
    It uses the data files: 
        -'IV_curve_y_0.1_k_0.01_G_0.001_T_0_r_0.0(125)_steps_120_tmax_40.npy'
        -'TUR_withQ_y_0.1_k_0.01_G_0.001_T_0.0005_steps_500_tmax_5000_nsamp_80000.npy'
    
    