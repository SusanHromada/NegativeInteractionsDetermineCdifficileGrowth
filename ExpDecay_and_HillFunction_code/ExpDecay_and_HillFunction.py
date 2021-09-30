#Written by Susan Hromada, Venturelli lab 2021

import numpy as np
from scipy.optimize import curve_fit

#Exp decay
def func(x, a, b):
    return a * np.exp(-b * x) 
popt, pcov = curve_fit(func, richness_data, OD_data, bounds=([0.01,0], [5,1]), maxfev=2000)

#Hill
def hill(x, a, n, EC50):
    return ((x**n)/((x**n)+(EC50**n)))
def return_hill_fit(xdata, ydata, p0):
    popt, pcov = curve_fit(hill, xdata, ydata, p0, bounds=([0,1,0],[100,10,100]),max_nfev=100000)
    return popt,pcov
popt, pcov = return_hill_fit(densities, OD_data_norm, [1,1,30])

