#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 13:01:57 2016

@author: bartonjo

This is a function to return value of the integral complementary error 
function to a power, n=1, 2, or 3: ierfc() for integrating once, 
i2erfc() for integrating twice, and i3erfc() for integrating thrice.

INPUTS
    Arguments
        x: A floating number or array that will be the domain of the 
        function
    Keywords
        n: The number of times erfc(x) will be integrated.  If n is 
        left blank, then it defaults to n=1.

RETURNS
    y: The solution of inerfc(x) in array format, even if x is scalar.
        
EXAMPLE
>>> from int_erfc import inerfc
>>> y = inerfc(x, n=2)     # y is set to the double integral of erfc(x)

 
"""
import numpy as np
from scipy.special import erf, erfc
from scipy.integrate import quad, odeint

# Numerically calculate i^nerfc(x) to n = 1, 2, or 3
def inerfc(x, n=1):
    # Check to see if x is a scalar. If so, make it an array.
    single = 0
    if not isinstance(x, np.ndarray):
        x = np.array([x])
        single = 1
    
    if all(x == 0):
        y = x
        n = 0
    
    if n == 1:
        y = []
        for i in range(len(x)):
            z = x[i]
            result = quad(erfc, 0, z)
            y.append(result[0])

        y = np.array(y)
        
    elif n == 2:
        if single == 1:
            x = np.asscalar(x)
            t = np.arange(0., x, x/100.)
            
        else:
            t = x
            
        y0 = [0., 0.]                   # ICs for [f', f]
        def func(y, t):
            return [erfc(t), y[0]]      # y' vector = [f'', f']
        
        y = odeint(func, y0, t)            
        y = y[:,1]
        if single == 1:
            y = np.array([y[-1]])       # only want the last element
        
    elif n == 3:
        if single == 1:
            x = np.asscalar(x)
            t = np.arange(0., x, x/100.)
            
        else:
            t = x
            
        y0 = [0., 0., 0.]       # ICs for [f'', f', f]
        def func(y, t):
            return [erfc(t), y[0], y[1]]    # y' vector = [f'', f']
        
        y = odeint(func, y0, t)            
        y = y[:,2]
        if single == 1:
            y = np.array([y[-1]])       # only want the last element
    
    return y


# Analytical solution for ierfc(x)
def i1erfc(x):    
    y = x*erfc(x) - np.exp(-x**2)/np.sqrt(np.pi) #+ 1/np.sqrt(np.pi)
    # Check to see if x is a scalar. If so, make y an array.
    if not isinstance(x, np.ndarray):
        y = np.array([y])
    
    return y
    

# Analytical solution for i2erfc(x)
def i2erfc(x):    
    y = x**2/2.*erfc(x) - x/2./np.sqrt(np.pi)*np.exp(-x**2) - \
            erf(x)/4.# + x/np.sqrt(np.pi))
    # Check to see if x is a scalar. If so, make y an array.
    if not isinstance(x, np.ndarray):
        y = np.array([y])
    
    return y
    
    
# Analytical solution for i3erfc(x)
def i3erfc(x):
    y = x**3*erfc(x)/6. - \
            (x**2 + 1)*np.exp(-x**2)/6./np.sqrt(np.pi) - \
            x*erf(x)/4. #+ \
            #x**2/2/np.sqrt(np.pi) + \
            #1/6./np.sqrt(np.pi)
    # Check to see if x is a scalar. If so, make y an array.
    if not isinstance(x, np.ndarray):
        y = np.array([y])
    
    return y    
    
    
    
    
    
    
    