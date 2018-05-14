# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:13:56 2017

@author: bartonjo

Compare TC and LP heat flux data
"""
import numpy as np
import matplotlib.pyplot as plt
import import_hfdecon_data as hfdata

# data attributes from files
shot = 167353#167476
tcname = 9#10
pname = 51#19

s = hfdata.getdecondata(shot=shot,tcname=tcname,pname=pname)
#----------------------------------------------------------------------
# data manipulation and plotting below
#

# calculate angle when heat flux is constant
tmin = 2500.
tmax = 4000.
zmin = np.min(np.where(s.lpt>tmin))
zmax = np.min(np.where(s.lpt>tmax))
av_ang = np.mean(s.ang[zmin:zmax])
# calculate tc heat flux using the finite slab approximation
    # Material constants for Mo
rho = 10.2e3                     # density [kg/m3]
cp = .23e3                       # specific heat [J/kg/K]
k = 138.                         # thermal conductivity [W/m/K]
kappa = k/rho/cp                # thermal diffusivity [m2/s]
d = .75e-2                      # measurement depth [m]
tcq = s.tcderiv*rho*cp*d

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(s.lptmed/1e3+.4,s.lpqmed/1e2/np.sin(np.radians(av_ang)), 
         'ks', markersize=10, fillstyle='none', 
         label='LP heat flux shifted 400 ms')
ax1.plot(s.tct[0::100]/1e3,tcq[0::100]/1e6/np.sin(np.radians(av_ang)), 
         'k.', markersize=15, label='TC heat flux')
ax2.plot(s.tct[0::100]/1e3,s.tcTemp[0::100], 'b', linewidth=2)

ax1.set_ylim([0,150])
ax2.set_ylim([100,700])

ax1.set_xlabel('Time [s]', fontsize=18)
ax1.set_ylabel(r'q$_{\mathrm{||}}$ [MW/m$^2$]', fontsize=18)
ax2.set_ylabel('Temperature [C]', color='b', fontsize=18)
ax2.tick_params(axis='y',colors='b')

h = '{:0.1f}'.format(av_ang)
plt.title('Heat flux data for shot 167476\n' +\
        r'P$_{\mathrm{tot}}$ = 10MW   $\theta$ = ' + h +\
        '$^o$ (angle to LP)')
ax1.legend(loc='upper left')
plt.show()