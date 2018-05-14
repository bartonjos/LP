# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:49:57 2017

@author: bartonjo

Calculate the ratio of surface temperature to end temperature using 
the solution of the finite slab 1-D heat equation for various values
of time.
"""
import numpy as np
from scipy.special import erfc

#---------------------------------------------------------------------
# Define parameters
#

# Material constants for pure Tungsten
rho = 19.25e3                   # density [kg/m3]
cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
k = 173.                        # thermal conductivity [W/m/K]
kappa = k/rho/cp                # thermal diffusivity [m2/s]


# User defined parameters and array initializations
a = 1.9e-3                  # slab thickness[m] 1.5flat 1.9mid 3.2top
t = range(1,7)              # seconds
t = np.array(t)*.5
num = range(100)             # number of terms in summation

#---------------------------------------------------------------------
# Use the solution of linear heat eaquation with constant flux on 
# surface and zero flux at x = a to calculate the ratio of Tsurf to 
# T(x=a), the bottom of the slab
#

# define ierfc function, which is the integral of erfc
def ierfc(z):
    y = z * erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)
    return y

# eval the summation part of the solution at x=0 and x=a
s0, s1 = 0., 0.
for i in num:
    # summation for surface, s0
    arg1 = (2*a*(1+i))/(2*np.sqrt(kappa*t))
    arg2 = (a*i)/np.sqrt(kappa*t)
    s0 += ierfc(arg1) + ierfc(arg2)
    # summation for T at x=a (i.e. the bottom of the slab), s1
    arg = (2*a*i + a)/(2*np.sqrt(kappa*t))
    s1 += 2*ierfc(arg)

# solution when x=a is defined as the surface
#for i in num:
#    # summation for surface, s1
#    arg1 = ((2*i+1)*a-a)/(2*np.sqrt(kappa*t))
#    arg2 = ((2*i+1)*a+a)/(2*np.sqrt(kappa*t))
#    s1 += ierfc(arg1) + ierfc(arg2)
#    # summation for T at x=0 (i.e. the bottom of the slab), s0
#    arg1 = ((2*i+1)*a)/(2*np.sqrt(kappa*t))
#    arg2 = ((2*i+1)*a)/(2*np.sqrt(kappa*t))
#    s0 += ierfc(arg1) + ierfc(arg2)
    

# dTsurf = dTbottom * R
R = s0/s1
print(R)