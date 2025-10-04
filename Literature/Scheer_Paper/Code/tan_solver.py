#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:41:42 2023

@author: david
"""

"""
This code solves the cubic equation (t^2+a)(t+1)=b given by the synchronization condition.
In addition it checks the stability of the obtained real solutions in order to obtain
 a full picture of the steady state.
"""

import numpy as np
import matplotlib.pyplot as plt 

def pv(a):
    return 1-a**2/3

def qv(a,b):
    return b-2*a*(1+a**2/9)/3

def isstable(t,b):
    return 2*b*t<=(1+t**2)**2#4*b*t+(1+t**2)**2>=0#

#Check stability for a given list of solutions
def classify(solutions,b):
    stable=[]
    unstable=[]
    for t in solutions:
        if isstable(t,b):
            stable.append(t)
        else:
            unstable.append(t)
    return stable, unstable

# Obtain all real solutions to the equation grouped by stability
def tan_solution(a,b):
    b=-b
    p=pv(a)
    q=qv(a,b)
    
    d=-(4*p**3+27*q**2) # calculate the discrminant
    
    if d==0: #three real solutions with one double root
        if p!=0:
            t1=3*q/p+a/3
            t2=-3/2*q/p+a/3
            return classify([t1,t2],b)
        else:
            t=a/3
            return list(t),list()
    elif d<0: #one real solution
        if p<0:
            t=a/3-2*np.sign(q)*np.sqrt(-p/3)*np.cosh(np.arccosh(-3*np.abs(q)/(2*p)*np.sqrt(-3/p))/3)
        elif p>0:
            t=a/3-2*np.sqrt(p/3)*np.sinh(np.arcsinh(3*q/(2*p)*np.sqrt(3/p))/3)
        else:
            t=a/3-q**(1/3)
        return classify([t],b)  
    else: #three real solutions
        solutions=[0,1,2]
        for k in range(3):
            solutions[k]=a/3+2*np.sqrt(-p/3)*np.cos(np.arccos(3*q/(2*p)*np.sqrt(-3/p))/3-2*np.pi*k/3)
        return classify(solutions,b)
 
# Obtain the solution as a function of a for a given b
def calculate_curves(alist,b):
    stable_1=[]
    avals_1=[]
    stable_2=[]
    avals_2=[]
    unstable_1=[]
    
    avals_u=[] # Give the range of , in which the different solutions are valid 
    crossed=False# Indicate whether left or right from the hysteretic region in case of only 1 solution
    decider=0
    for a in alist:
        stable,unstable=tan_solution(a,b)
        if len(unstable)==0 and not crossed:
            stable_1.append(stable[0])
            avals_1.append(a)
        elif len(unstable)==0 and crossed:
            stable_2.append(stable[0])
            avals_2.append(a)
        else:
            if len(stable)<=1:#this is a fix to identify the actually stable solution
                continue
                
            if not crossed:
                crossed=True
                try:    
                    decider=int(np.abs(stable[0]-stable_1[-1])>np.abs(stable[1]-stable_1[-1]))
                except:
                    print(stable_1)
            #this identification is in the wrong order but it works. Keep until problem with stability identification is fixed
            stable_1.append(stable[decider])
            avals_1.append(a)
            stable_2.append(stable[1^decider])
            avals_2.append(a)
            unstable_1.append(unstable[0])
            avals_u.append(a)
    return [np.array(stable_1),np.array(stable_2),np.array(unstable_1)],[np.array(avals_1),np.array(avals_2),np.array(avals_u)]
