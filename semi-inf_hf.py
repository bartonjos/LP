# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:52:57 2017

@author: bartonjo

Calculate the surface temperature with constant and transient heat 
fluxes on the surface and semi-infinite BCs
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

#---------------------------------------------------------------------
# Define parameters
#

# Material constants for pure Tungsten samples
rho = 19.25e3                   # density [kg/m3]
cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
k = 173.                        # thermal conductivity [W/m/K]
kappa = k/rho/cp                # thermal diffusivity [m2/s]
a = 1.5e-3                      # flat sample thickness [m]
am = 1.9e-3                     # mid way up angled thickness [m]
at = 3.2e-3                     # all way up angled thickness [m]

# User defined parameters and array initializations
t = np.linspace(1e-4,3,10000)              # seconds
x = 0                               # [m] from surface
F = np.array([1e6,8e6,21e6])        # flux on surface [W/m2]
phip = np.array([2e6,16e6,33e6])      # transient flux power density
elmtime = np.array([5e-3,5e-3,5e-4])      # [s] ELM duration
elmt = np.zeros(np.size(elmtime))             # initialize array
for j in range(3):
    for i in t:
        if (i<=elmtime[j]):
            elmt[j] = i      # put the elm duration on t scale
phi = 2*phip*elmtime                # tran flux energy density [J/m2]
elmf = np.array([10,10,170])            # "ELM" frequency Hz
numelms = np.max(t)*elmf
# time when trans happens [s]
for j in range(3):
    if j<2:
        tau1 = [t[i] for i in range(np.size(t)) \
            if np.mod(i+1,np.round(np.size(t)/numelms[j]))==0]
    else:
        tau2 = [t[i] for i in range(np.size(t)) \
            if np.mod(i+1,np.round(np.size(t)/numelms[j]))==0]

Tempt = np.zeros((3, np.size(t)))   # temperature [K]
Telms = np.zeros((1, np.size(t)))
Telms2 = np.zeros((1, np.size(t)))
num = range(50)                     # number of terms in summation


#---------------------------------------------------------------------
# Evaluate solution of linear heat eaquation with constant flux on 
# surface and zero flux at x = inf
#

# calculate the shifted time for the elm heat flux
def tshift(t,s,telm):
    shift = [i for i in t if i>=telm]
    z = [np.min(t)/10. for i in t if (i-s)<=0]
    z.extend(shift)
    z = np.array(z[0:np.size(t)])
    return z

for j in range(3):
    if j<2:
        # eval the temp vs time with a constant background heat flux, F
        Tinter = F[j]/k*(2*np.sqrt(kappa*t/np.pi)*np.exp(-x**2/4/kappa/t) - 
                x*erfc(x/2/np.sqrt(kappa*t)))    
        # evaluate the temp vs time for ELMs
        for m in range(np.size(tau1)):
            ts = tshift(t,tau1[m],elmt[j])
            Telms += phi[j]/k*np.sqrt(kappa/np.pi/ts)*np.exp(-x**2/4/kappa/ts) * \
                    np.piecewise(ts,
                    [ts==np.min(t)/10., ts>np.min(t)/10.],[0,1])
        
        # Linearity allows us to add the inter-ELM and ELM solutions
        Tempt[j][:] = Tinter + Telms
    else:
        # eval the temp vs time with a constant background heat flux, F
        Tinter2 = F[j]/k*(2*np.sqrt(kappa*t/np.pi)*np.exp(-x**2/4/kappa/t) - 
                x*erfc(x/2/np.sqrt(kappa*t)))    
        # evaluate the temp vs time for ELMs
        for m in range(np.size(tau2)):
            ts2 = tshift(t,tau2[m],elmt[j])
            Telms2 += phi[j]/k*np.sqrt(kappa/np.pi/ts2)*np.exp(-x**2/4/kappa/ts2) * \
                    np.piecewise(ts2,
                    [ts2==np.min(t)/10., ts2>np.min(t)/10.],[0,1])
        
        # Linearity allows us to add the inter-ELM and ELM solutions
        Tempt[j][:] = Tinter2 + Telms2

#---------------------------------------------------------------------
# Calculate the slope of the temperature on the surface from the slope
# of the temperature at a distance x=a. How much greater is the 
# slope on the surface?
if x==0:
    x=a
t_shot = range(1,7)
t_shot = np.array(t_shot)*.5
ms = np.sqrt(4*kappa*t_shot/np.pi)
ml = np.sqrt(4*kappa*t_shot/np.pi)*np.exp(-x**2/4/kappa/t_shot) - \
        x*erfc(x/2/np.sqrt(kappa*t_shot))
R = ms/ml
#print(R)
        

#---------------------------------------------------------------------
# Plot solutions
#

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

#plt.plot(oldt, oldTemp, linewidth=3, label='flush')
plt.plot(t+2, (Tempt[0][:]+Toffset), linewidth=3, 
         label='flush,'+ 
         ' ref 167297: Pnbi=3.1MW  ELMfreq=10Hz')
h = '{:0.0f}'.format(15)
plt.plot(t+2, (Tempt[1][:]+Toffset), linewidth=3, 
         label=h + ' degree angle,' + 
         ' ref 167297: Pnbi=3.1MW  ELMfreq=10Hz')
h = '{:0.0f}'.format(15)
plt.plot(t+2, (Tempt[2][:]+Toffset), linewidth=3, 
         label=h + ' degree angle,' +
         ' ref 167318: Pnbi=7.3MW  ELMfreq=170Hz')

plt.plot([0, np.max(t)+2], [3422, 3422], 'k')
plt.plot([0, np.max(t)+2], [1000, 1000], 'k')
plt.plot([0, np.max(t)+2], [0,0], 'k')

#plt.xlim([0, np.floor(np.max(t))*1.05])
#plt.ylim([np.min(Tempt)*.95, 3695*1.15])
plt.xlim([1.2,5.2])
plt.ylim([50, 3700])

#plt.title('Constant ' + h + ' MW/m$^2$ of applied heat flux\n' + 
#        'Transient ' + h1 + ' MW/m$^2$ heat flux at ' + h2 + ' Hz'
#        , fontsize=18)
#plt.title('158463 (He), Ptot = 2.7MW, Pech = 2.3MW\n' + 
#        'Ip = 1.09, Bt = -1.97'
#        , fontsize=18)
#plt.title('166718 (D), Ptot = 4.6MW, Pnbi = 4.5MW\n' + 
#        'Ip = 1.13, Bt = -2.07, ELMfreq = ~35Hz'
#        , fontsize=18)
#plt.title('167297, Ptot = 3.6MW, Pnbi = 3.1MW\n' + 
#        'Ip = 1.29, Bt = -2.06, ELMfreq = ~10Hz'
#        , fontsize=18)
plt.xlabel('Time [s]', fontsize=18)
plt.tick_params(axis='x', labelsize=16)
plt.ylabel('Surface temperature [C]', fontsize=18)
#plt.ylabel('T-T$_o$ [K]', fontsize=18)
plt.tick_params(axis='y', labelsize=16)
plt.legend(loc='upper left')
fig.subplots_adjust(top=0.89)
plt.show()