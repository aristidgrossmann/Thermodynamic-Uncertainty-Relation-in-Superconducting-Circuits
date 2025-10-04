import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import LogLocator

import os
import sys

mpl.rcParams.update({
    "text.usetex": False,            # Use matplotlib's internal math rendering
    "mathtext.fontset": "cm",        # Use Computer Modern (LaTeX-like)
    "font.family": "serif",          # Serif font family
    "font.size": 11,                 # Base font size
    "xtick.labelsize": 10,           # Tick label font size
    "ytick.labelsize": 10,
    "axes.labelsize": 11,            # Axis label font size
    "legend.fontsize": 10,           # Legend font size (default: 10)
})


os.chdir('C:/Users/arist/Desktop/SS25/BA/Abgabe/figures/04/JJ')

sys.path.append( os.getcwd())

from Exact_Solutions import  I_tilde_minus_plus, I_tilde_plus_minus, T_1_paper






i_0 = 0.8
N_phi = 1000
phi_vals = np.linspace(-np.pi/2, 1.5*np.pi, N_phi)

D_vals = np.array([0.1, 0.5, 1.5])


for i in range(len(D_vals)):
    PDF_minus = np.ones_like(phi_vals)
    PDF_plus = np.ones_like(phi_vals)

    MFPT = T_1_paper(i_0 = i_0, D = D_vals[i])*D_vals[i]*(1-np.exp(-i_0*2*np.pi/D_vals[i]))

    for j in range(len(PDF_minus)): 
        PDF_minus[j] = I_tilde_minus_plus(phi_vals[j], i_0 = i_0, D = D_vals[i])/MFPT
        PDF_plus[j] = I_tilde_plus_minus(phi_vals[j], i_0 = i_0, D = D_vals[i])/MFPT

        if j % 10 == 0:
            print(j)

    filename = "PDF_minus_" + str(D_vals[i]) + ".txt"
    np.savetxt(filename, PDF_minus)

    filename = "PDF_plus_" + str(D_vals[i]) + ".txt"
    np.savetxt(filename, PDF_plus)









