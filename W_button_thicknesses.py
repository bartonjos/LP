#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code will help in the design of the W button thicknesses to 
calculate how hot they will get with a given heat flux.

Created on Wed Jan 11 13:16:07 2017

@author: bartonjo
"""

import numpy as np
import matplotlib.pyplot as plt
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
a = np.linspace(1e-4,3e-3)          # slab thickness [m]
a1 = 1.5e-3
a2 = 2e-3
a3 = 3e-3
t = np.linspace(1e-3,3,10000)              # seconds
t1 = 0.5
t2 = 1.
t3 = 2.
x = 0                               # [m] from surface
x1 = np.linspace(0,a2)
F = 8e6                             # flux on surface [W/m2]
phip = 16e6                           # transient flux power density
elmtime = 5e-3                       # [s] ELM duration
for i in t:
    if (i<=elmtime):
        elmt = i                    # put the elm duration on t scale
phi = 2*phip*elmt                   # tran flux energy density [J/m2]
elmf = 10                           # "ELM" frequency Hz
numelms = np.max(t)*elmf
# time when trans happens [s]
tau = [t[i] for i in range(np.size(t)) \
        if np.mod(i+1,np.round(np.size(t)/numelms))==0]

Tempa = np.zeros((3, np.size(a)))   # temperautre [K]
Tempt = np.zeros((3, np.size(t)))
Telms = np.zeros((1, np.size(t)))
num = range(10)                     # number of terms in summation

#---------------------------------------------------------------------
# Evaluate solution of linear heat eaquation with constant flux on 
# surface and zero flux at x = a
#

# define ierfc function, which is the integral of erfc
def ierfc(z):
    y = z * erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)
    return y

# define the heavyside function for transient heat fluxes
def heavy(t,tau,elmt):
    H = np.zeros(np.size(t))
    for s in tau:
        p1 = np.piecewise(t, [(t-s)<=0, (t-s)>0], [0, 1])
        p2 = np.piecewise(t, [(t-s-elmt)<=0, (t-s-elmt)>0], [0, 1])
        H += p1-p2
    return H

# convert elm pulse over dt to instantaneous pulse
def convertphi(phi, s, dt):
    P = phi/s
#    print('{:0.2e}'.format(np.sqrt(dt*tmin)))
    return P

# calculate the shifted time for the elm heat flux
def tshift(t,s,elmt):
    shift = [i for i in t if i>=elmt]
    z = [np.min(t)/10. for i in t if (i-s)<=0]
    z.extend(shift)
    z = np.array(z[0:np.size(t)])
    return z

H = heavy(t,tau,elmt)
j = 0
for l in [a1, a1]:#, a3]: 
    # eval the temp vs time for each thickness value
#    s = 0.
#    for i in num:
#        arg1 = (2*l*(1+i) - x)/(2*np.sqrt(kappa*t))
#        arg2 = (2*l*i + x)/(2*np.sqrt(kappa*t))
#        s += ierfc(arg1) + ierfc(arg2)
#    Tinter = s * (-2)*np.sqrt(kappa*t)/k * F
    Tinter = F/k*(2*np.sqrt(kappa*t/np.pi)*np.exp(-x**2/4/kappa/t) - 
            x*erfc(x/2/np.sqrt(kappa*t)))    
    
    # evaluate the temp vs time for ELMs
    # first convert elm fluxes
    s1 = 0.
    for i in num:
        arg1 = -(2*l*(i+1))**2/4./kappa/elmt
        arg2 = -(2*l*i)**2/4./kappa/elmt
        s1 += np.exp(arg1) + np.exp(arg2)
    P = convertphi(phi,s1,elmt)
#    P = phi
#    print('{:.0f}'.format(P))
    if (j==0):
        P=0
    
    # Now evaluate the temp vs time for ELMs
    for m in range(np.size(tau)):
        ts = tshift(t,tau[m],elmt)
#            ts = t[n] - tau[m]
        s1 = 0.
        for i in num:
            arg1 = -(2*l*(i+1) - x)**2/4./kappa/ts
            arg2 = -(2*l*i + x)**2/4./kappa/ts
            s1 += np.exp(arg1) + np.exp(arg2)
#        s1 = np.exp(-x**2/4/kappa/ts)
        Telms += P/np.sqrt(rho*cp*k*np.pi*ts) * s1 * \
                np.piecewise(ts,
                [ts==np.min(t)/10., ts>np.min(t)/10.],[0,1])
    
    # Linearity allows us to add the inter-ELM and ELM solutions
    Tempt[j][:] = Tinter + Telms
    j += 1

#### iterate over time
#j = 0
#for l in [a1]: 
#    # evaluate temp at each time step
#    for n in range(1, np.size(t)):
#        s, Treg, Telms= 0., 0., 0.
#        # eval the temp from inter-ELM heat flux
#        for i in num:
#            arg1 = (2*l*(1+i) - x)/(2*np.sqrt(kappa*t[n]))
#            arg2 = (2*l*i + x)/(2*np.sqrt(kappa*t[n]))
#            s += ierfc(arg1) + ierfc(arg2)
#        Treg = s * (-2)*np.sqrt(kappa*t[n])/k * F
#        # evaluate the temp from ELMs
#        for m in range(np.size(tau)):
#            ts = t[n] - tau[m]
#            if (ts>0):
#                s1 = 0.
#                for i in num:
#                    arg1 = -(2*l*(i+1) - x)**2/4./kappa/ts
#                    arg2 = -(2*l*i + x)**2/4./kappa/ts
#                    s1 += np.exp(arg1) + np.exp(arg2)            
#                Telms += phi/np.sqrt(rho*cp*k*np.pi*ts) * s1
#        Tempt[j][n] = Treg + Telms + Tempt[j][n-1]
#    j += 1
####

#j = 0
#for p in [t1,t2,t3]: 
#    # eval the temp vs thickness for each time value
#    s = 0.
#    for i in num:
#        arg1 = (2*a*(1+i) - x)/(2*np.sqrt(kappa*p))
#        arg2 = (2*a*i + x)/(2*np.sqrt(kappa*p))
#        s += ierfc(arg1) + ierfc(arg2)
#    Tempa[j][:] = -2*F*np.sqrt(kappa*p)/k * s
#    j += 1
    
#j = 0
#for p in [t1,t2,t3]: 
#    # eval the temp vs depth for each time value
#    s = 0.
##    for n in num:
##        i = n + 1
##        arg = (-1)**i/i**2 * np.exp(-kappa*i**2*np.pi**2*p/a2**2) * \
##                np.cos(i*np.pi*x1/a2)
##        s += arg
##    Tempa[j][:] = F*p/rho/cp/a2 + F*a2/k * \
##                ((3*x1**2 - a2**2)/(6*a2**2) - 2/np.pi**2*s)
#    for i in num:
#        arg1 = (2*a2*(1+i) - x1)/(2*np.sqrt(kappa*p))
#        arg2 = (2*a2*i + x1)/(2*np.sqrt(kappa*p))
#        s += ierfc(arg1) + ierfc(arg2)
#    Tempa[j][:] = -2*F*np.sqrt(kappa*p)/k * s
#    j += 1


#---------------------------------------------------------------------
# Plot solutions
#

# Temperature vs thickness for 3 specific times
#fig, ax = plt.subplots()
#h = '{:0.1f}'.format(t1)
#plt.plot(a*1e3, Tempa[0][:]-273, linewidth=3, label=h + ' s')
#h = '{:0.1f}'.format(t2)
#plt.plot(a*1e3, Tempa[1][:]-273, linewidth=3, label=h + ' s')
#h = '{:0.1f}'.format(t3)
#plt.plot(a*1e3, Tempa[2][:]-273, linewidth=3, label=h + ' s')
#
#plt.xlim([np.min(a)*1e3*.8, np.max(a)*1e3*1.2])
#plt.ylim([20, np.floor(np.max(Tempa)-273)*1.2])
#
#h = '{:0.1f}'.format(F/1e6)
#plt.title('Constant ' + h + ' MW/m$^2$ of applied heat flux',
#          fontsize=18)
#plt.xlabel('Slab thickness [mm]', fontsize=18)
#plt.tick_params(axis='x', labelsize=16)
#plt.ylabel('Temperature [C]', fontsize=18)
#plt.tick_params(axis='y', labelsize=16)
#plt.legend(loc='best')
#plt.show()

# Temperature vs depth for 3 specific times
#fig, ax = plt.subplots()
#h = '{:0.1e}'.format(t1)
#plt.plot(x1*1e3, Tempa[0][:], linewidth=3, label=h + ' s')
#h = '{:0.1e}'.format(t2)
#plt.plot(x1*1e3, Tempa[1][:], linewidth=3, label=h + ' s')
#h = '{:0.1e}'.format(t3)
#plt.plot(x1*1e3, Tempa[2][:], linewidth=3, label=h + ' s')
#
#plt.xlim([0, np.max(x1)*1e3*1.05])
##plt.ylim([np.floor(np.min(Tempa)-273)*1.05, 
##          np.floor(np.max(Tempa)-273)*1.05])
#
#h = '{:0.1f}'.format(F/1e6)
#plt.title('Constant ' + h + ' MW/m$^2$ of applied heat flux',
#          fontsize=18)
#plt.xlabel('Depth [mm]', fontsize=18)
#plt.tick_params(axis='x', labelsize=16)
#plt.ylabel('Temperature [$\Delta$T]', fontsize=18)
#plt.tick_params(axis='y', labelsize=16)
#plt.legend(loc='best')
#plt.show()


# Temperature vs time for 3 specific thicknesses
Toffset = 100-273+293
fig, ax = plt.subplots()
#h = '{:0.1f}'.format(a1*1e3)
#plt.plot(t+2, Tempt[0][:]+Toffset, linewidth=3, label=h + ' mm')
#h = '{:0.1f}'.format(a1*1e3)
#plt.plot(t+2, Tempt[1][:]+Toffset, linewidth=3, 
#label=h + ' mm with ELMs')
#h = '{:0.1f}'.format(a3*1e3)
#plt.plot(t, Tempt[2][:], linewidth=3, label=h + ' mm')

plt.plot(oldt, oldTemp, linewidth=3, label='flush')
h = '{:0.0f}'.format(15)
plt.plot(t+2, Tempt[1][:]+Toffset, linewidth=3, 
         label=h + ' degree angle')

plt.plot([0, np.max(t)+2], [3422, 3422], 'k')
plt.plot([0, np.max(t)+2], [1000, 1000], 'k')
plt.plot([0, np.max(t)+2], [0,0], 'k')

#plt.xlim([0, np.floor(np.max(t))*1.05])
plt.ylim([np.min(Tempt)*.95, 3695*1.15])
plt.xlim([1.2,5.2])
#plt.ylim([50, 500])

h = '{:0.1f}'.format(F/1e6)
h1 = '{:0.1f}'.format(phip/1e6)
h2 = '{:0.0f}'.format(elmf)
#plt.title('Constant ' + h + ' MW/m$^2$ of applied heat flux\n' + 
#        'Transient ' + h1 + ' MW/m$^2$ heat flux at ' + h2 + ' Hz'
#        , fontsize=18)
#plt.title('158463 (He), Ptot = 2.7MW, Pech = 2.3MW\n' + 
#        'Ip = 1.09, Bt = -1.97'
#        , fontsize=18)
#plt.title('166718 (D), Ptot = 4.6MW, Pnbi = 4.5MW\n' + 
#        'Ip = 1.13, Bt = -2.07, ELMfreq = ~35Hz'
#        , fontsize=18)
plt.title('167297, Ptot = 3.6MW, Pnbi = 3.1MW\n' + 
        'Ip = 1.29, Bt = -2.06, ELMfreq = ~10Hz'
        , fontsize=18)
plt.xlabel('Time [s]', fontsize=18)
plt.tick_params(axis='x', labelsize=16)
plt.ylabel('Surface temperature [C]', fontsize=18)
#plt.ylabel('T-T$_o$ [K]', fontsize=18)
plt.tick_params(axis='y', labelsize=16)
plt.legend(loc='upper left')
fig.subplots_adjust(top=0.89)
plt.show()
