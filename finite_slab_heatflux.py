#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code is used to calculate the conditions necessary for the MTX 
tungsten coated moly inserts to bow and melt.

Created on Tue Sep  6 15:30:47 2016

@author: bartonjo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from int_erfc import inerfc

#---------------------------------------------------------------------
# Define parameters
#

# Material constants for Moly inserts
rho = 10.2e3                    # density [kg/m3]
cp = 0.23e3                     # specific heat [J/kg/K]  0.23 J/g/K
a = 4.7e-3              # slab thickness [m] cham@4.7e-3 full@9.525e-3
k = 138.                        # thermal conductivity [W/m/K]
kappa = k/rho/cp                # thermal diffusivity [m2/s]

# User defined parameters
t = 3.0
t2 = np.linspace(10./1000.,3)                         # time [s]
x = a#7.5e-3            # distance from surface [m] probe @7.5e-3
F = np.linspace(0,80e6)         # heat flux to surface [W/m2]

Temp = np.zeros((2, np.size(F)))        # 2 depths and a flux range
Temp2 = np.zeros((2, np.size(t2)))      # 2 depths and a time range

num = range(10)                   # number of terms in summation

#---------------------------------------------------------------------
# Evaluate solution of linear heat eaquation with constant flux on 
# surface and zero flux at x = a
#

# define ierfc function, which is the integral of erfc
def ierfc(z):
    y = z * erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)
    return y

j = 0
for x1 in [0, x]:
    s = 0.0
    for i in num:
        arg1 = (2*a*(1+i) - x1)/(2*np.sqrt(kappa*t))
        arg2 = (2*a*i + x1)/(2*np.sqrt(kappa*t))
        s += ierfc(arg1) + ierfc(arg2)
#        print('{0:.2e}'.format(s) + '\n')
    Temp[j][:] = -2*F*np.sqrt(kappa*t)/k * s
    j += 1

#oldT = np.zeros((2, np.size(F)))
#j = 0
#for x1 in [0, a]:
#    s = 0.0
#    for i in num:
#        coeff = 2*i + 1
#        arg1 = (coeff * a - x1)/(2*np.sqrt(kappa*t))
#        arg2 = (coeff * a + x1)/(2*np.sqrt(kappa*t))
#        s += ierfc(arg1) + ierfc(arg2)
#    oldT[j][:] = -2*F/k*np.sqrt(kappa*t) * s
#    j += 1

# Evaluate the Temperature when the heat flux increases with time: i.e.
# F = f*t, where f is some constant
f = 4e6/min(t2)
j = 0
for x1 in [0, x]:
    s = 0.
    for i in num:
        arg1 = (2*a*(i+1) - x1)/2./np.sqrt(kappa*t2)
        arg2 = (2*a*i + x1)/2./np.sqrt(kappa*t2)
        s += inerfc(arg1, n=3) + inerfc(arg2, n=3)
        
    Temp2[j][:] = -8*(f*t2)*np.sqrt(kappa*t2)/k * s
    j += 1

#---------------------------------------------------------------------
# Plot solutions
#

# Temperature at x vs heat flux
fig = plt.subplots()
plt.plot(F/1e6, Temp[0][:]-273, linewidth=3, label='Surface')
h = '{0:.2f}'.format(x*100)
plt.plot(F/1e6, Temp[1][:]-273, linewidth=3, label='Depth = ' + h + \
         ' cm')

#plt.plot(F/1e6, oldT[1][:]-273, linewidth=3, label='Surface*')
#plt.plot(F/1e6, oldT[0][:]-273, linewidth=3, label='Depth* = ' + \
#       h + ' cm')

flat = np.array([0, max(F)])/1e6
plt.plot(flat, np.ones(2)*2623, label='Mo melts')
plt.plot(flat, np.ones(2)*3422, label='W melts')
plt.plot(flat, np.ones(2)*3457, label='Graphite vaporizes')

plt.xlim([min(F)/1e6, 30])
plt.ylim([0, 5e3])

h = '{0:0.1f}'.format(t)
plt.title('Temperature after ' + h + ' s of applied heat flux', \
            fontsize=18)
plt.xlabel('Heat flux deposited onto surface [MW/m$^2$]', fontsize=18)
plt.tick_params(axis='x', labelsize=16)
plt.ylabel('Temperature [C]', fontsize=18)
plt.tick_params(axis='y', labelsize=16)
plt.legend(loc='best')
folder = '/Users/bartonjo/Documents/Logs/metal rings campaign 2016/'
#plt.savefig(folder + 'T_vs_HF_at_edge.png')
plt.show()

# heat flux vs angle of incidence
fig = plt.subplots()
ang = np.linspace(0,15)        # angle in degrees
qpar = np.array([73., 122.])    # parallel heat flux
h1 = '{0:0}'.format(int(qpar[0]))
h2 = '{0:0}'.format(int(qpar[1]))
plt.plot(ang, qpar[0]*np.sin(np.radians(ang)), linewidth=3, \
        label=r'$q_{\parallel }$ = ' + h1 + ' MW/m$^2$')
plt.plot(ang, qpar[1]*np.sin(np.radians(ang)), linewidth=3, \
        label=r'$q_{\parallel }$ = ' + h2 + ' MW/m$^2$')
#plt.plot(np.array([min(ang), max(ang)]), np.array([qpar, qpar]), \
#            label='Parallel to B')

plt.xlim([0, 15])

#plt.title('Parallel heat flux = ' + '{0:0}'.format(int(qpar)) + \
#            ' [MW/m$^2$]', fontsize=18)
plt.xlabel(r'$ \theta $ [degrees]', fontsize=18)
plt.ylabel('Deposited heat flux [MW/m$^2$]', fontsize=18)
plt.tick_params(axis='x', labelsize=16)
plt.tick_params(axis='y', labelsize=16)
plt.legend(loc='best')
folder = '/Users/bartonjo/Documents/Logs/metal rings campaign 2016/'
#plt.savefig(folder + 'HF_vs_angle_' + '{0:0}'.format(int(qpar)) + \
#       '.png')
plt.show()

## Temperature vs time with applied constant heat flux
#fig = plt.subplots()
#plt.plot(t*1000, Temp[0][:]-273, linewidth=3, label='Surface')
#h = '{0:.2f}'.format(x*100)
#plt.plot(t*1000, Temp[1][:]-273, linewidth=3, label='Depth = ' + \
#       h + ' cm')
#
#h = '{0:0}'.format(int(F/1e6))
#plt.title('Temperature with ' + h + ' MW/m$^2$ of applied heat flux',\
#            fontsize=18)
#plt.xlabel('time [ms]', fontsize=18)
#plt.tick_params(axis='x', labelsize=16)
#plt.ylabel('Temperature [C]', fontsize=18)
#plt.tick_params(axis='y', labelsize=16)
#plt.legend(loc='best')
#plt.show()











