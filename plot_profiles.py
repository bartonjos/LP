# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 11:31:15 2018

@author: bartonjo

Import LP data from median filter profile dat files, plot them, fit
a curve to them, and save the curve
"""

import sys
s = '/Users/bartonjo/PyFiles/LP'
if s not in sys.path:
    sys.path.insert(0, s)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cookbook import savitzky_golay as smooth

# Get data from file into dataframe
jsat = pd.read_table('171433_jsat_med.dat', sep='\s+', header=None,
                     names=['R','jsat'])
qperp = pd.read_table('171433_qdata_med.dat', sep='\s+', header=None,
                     names=['R','qperp'])

# sort the data into ascending order in R and average redundant values
# jsat array sorting:
jsat.sort_values(by=['R'], inplace=True)
jsat_sort = np.array(jsat.jsat)
r_j = np.array(jsat.R)
rdum = np.array([])  # set up new array for R and jsat
jdum = np.array([])
for i in range(r_j.size):
    # find redundant R values
    x = r_j[i]
    if x not in rdum:
        c = 0  # index for number of the same R values
        for j in range(r_j.size):
            if r_j[j] == x:
                c += 1
        # save value of R and save average jsat at that R
        rdum = np.append(rdum, [x])
        jdum = np.append(jdum, [np.mean(jsat_sort[i+c-1])])
# qperp array sorting:
qperp.sort_values(by=['R'], inplace=True)
qperp_sort = np.array(qperp.qperp)
r_q = np.array(qperp.R)
rqdum = np.array([])  # set up new array for R and jsat
qdum = np.array([])
for i in range(r_q.size):
    # find redundant R values
    x = r_q[i]
    if x not in rqdum:
        c = 0  # index for number of the same R values
        for j in range(r_q.size):
            if r_q[j] == x:
                c += 1
        # save value of R and save average jsat at that R
        rqdum = np.append(rqdum, [x])
        qdum = np.append(qdum, [np.mean(qperp_sort[i+c-1])])

# convert to better units
rdum = rdum * 100  # m to cm
rqdum = rqdum * 100  # m to cm
e = 1.602e-19
jdum = jdum/e*1e4  # A/cm2 to ions/m2/s

# fit a smoothed curve to the data
jsat_smooth = smooth(jdum, 11, 3)
qperp_smooth = smooth(qdum, 11, 3)

# plot everything
plt.figure()
plt.subplot(2,1,1)
fsize = 22  # font size for axes
#plt.plot(rdum, jdum/1e22, 's', label='data')
plt.plot(jsat.R*100, jsat.jsat/e*1e4/1e22, 's', label='data')
plt.plot(rdum, jsat_smooth/1e22, label='fit')
    # sample locations relative to steady OSP at 146.7 cm
osploc = 146.7  # cm
plt.plot([147.376-osploc, 147.376-osploc+.6],[0.1,0.1], 'k')
plt.plot([148.247-osploc, 148.247-osploc+.6],[0.1,0.1], 'k')
plt.plot([149.103-osploc, 149.103-osploc+.6],[0.1,0.1], 'k')
plt.plot([0,0],[0,3],'k')
plt.xlim([-3,10])
plt.ylim([0,2])
plt.ylabel('particle flux\n' + r'[$\times 10^{22}$ ions/m$^2$/s]',
                                  fontsize=fsize)
plt.tick_params(labelsize=fsize-2)
plt.legend()

plt.subplot(2,1,2)
#plt.plot(rqdum, qdum, 's', label='data')
plt.plot(qperp.R*100, qperp.qperp, 's', label='data')
plt.plot(rqdum, qperp_smooth, label='fit')
    # sample locations relative to steady OSP at 146.7 cm
plt.plot([147.376-osploc, 147.376-osploc+.6],[0.03,0.03], 'k')
plt.plot([148.247-osploc, 148.247-osploc+.6],[0.03,0.03], 'k')
plt.plot([149.103-osploc, 149.103-osploc+.6],[0.03,0.03], 'k')
plt.plot([0,0],[0,1],'k')
plt.xlim([-3,10])
plt.ylim([0,.6])
plt.ylabel('heat flux\n' + r'[MW/m$^2$]', fontsize=fsize)
plt.xlabel('distance from OSP [cm]', fontsize=fsize)
plt.tick_params(labelsize=fsize-2)
plt.legend()
plt.subplots_adjust(left=.17, bottom=.13)

# save curve to file