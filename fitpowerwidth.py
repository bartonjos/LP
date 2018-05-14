#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fit SOL power width

This file will read in heat flux data profiles from the Langmiur probes and
fit an exponential to the right hand side of the profile, giving the
lambda_q (i.e. the SOL power width).  The script assumes the the filename
has the format ######_qdata.dat, so that the title can be the shot # and 
the file extention can be read by numpy.loadtxt. 

Example:
>>> run fitpowerwidth.py 167229_qdata.dat


@author: bartonjo
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#-----------------------------------------------------------------------------
# Make sure a file is being read; otherwise, return syntax message.
#
if (len(sys.argv) >1):
	infile = sys.argv[1]
 
else:
    sys.exit("syntax: run fitpowerwidth.py <infile>")
         
#-----------------------------------------------------------------------------
# Load the data into arrays
#
x, y = np.loadtxt(infile, unpack=True, usecols=(0, 1))

#-----------------------------------------------------------------------------
# Set up the plot of the raw x, y data
#
fig = plt.figure()
ax = fig.add_subplot(111)   # Placeholder in case you want additional plots
ax.plot(x, y, 'r.', ms=10., label='Probe Data')
xrange = max(x) - 0.
yrange = 200. - min(y)
minx = 0. - 0.05*xrange
maxx = max(x) + 0.05*xrange
miny = min(y) - 0.05*yrange
maxy = 200. + 0.05*yrange
ax.axis([minx, maxx, miny, maxy])
ax.tick_params(axis='x', labelsize=24)
ax.tick_params(axis='y', labelsize=24)
plt.title('Shot ' + infile[:6], fontsize=28)
plt.xlabel('$R-R_{sep}\ (m)$', fontsize=28)
plt.ylabel('$q_{\parallel}\ (MW/m^2)$ ', fontsize=28)

#-----------------------------------------------------------------------------
# Find the x and y values that will be used for the fitting
#
win = np.where(x > 0.)
xwin = x[win]
ywin = y[win]
#plt.plot(xwin,ywin,'bs')

#-----------------------------------------------------------------------------
# Do the fitting calculations
#
func = lambda z, a, b, c: a + b*z * np.exp(-c*z)    # function to be fitted

sigma = np.ones(len(xwin))  # weights of each point to be fit
p0 = np.ones(3)*10          # initial guesses for a, b, c

p, pcov = curve_fit(func,xwin,ywin,p0,sigma)    # p: coeffs, pcov: covar matrx

perr = np.sqrt(np.diag(pcov))                   # Standard Errors calculation

# R^2 calculation:
residuals = ywin - func(xwin, p[0], p[1], p[2])
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ywin - np.mean(ywin))**2)
r_sq = 1. - (ss_res/ss_tot)

xfit = np.linspace(0., 0.2, 100)
yfit = func(xfit, p[0], p[1], p[2])

print('\nBest fit:')
print(u'a = {:g} \u00B1 {:g}'.format(p[0],perr[0]))
print(u'b = {:g} \u00B1 {:g}'.format(p[1],perr[1]))
print(u'c = {:g} \u00B1 {:g}'.format(p[2],perr[2]))

#-----------------------------------------------------------------------------
# Now compare to a pure exponential fit starting from the max of the prev fit
#
fitwin = np.where(yfit == max(yfit))
win2 = np.where(x >= xfit[fitwin[0]])
xwin2 = x[win2]
ywin2 = y[win2]

#func2 = lambda z, a, b, c: a + b * np.exp(-c*z)     # function to be fitted
func2 = lambda z, b, c: b * np.exp(-c*z)     # function to be fitted

q, qcov = curve_fit(func2, xwin2, ywin2)

qerr = np.sqrt(np.diag(qcov))               # Standard Errors calculation

# R^2 calculation:
#residuals2 = ywin2 - func2(xwin2, q[0], q[1], q[2])
residuals2 = ywin2 - func2(xwin2, q[0], q[1])
ss_res2 = np.sum(residuals2**2)
ss_tot2 = np.sum((ywin2 - np.mean(ywin2))**2)
r_sq2 = 1. - (ss_res2/ss_tot2)

xfit2 = np.linspace(0.,0.2,100)
#yfit2 = func2(xfit2,q[0],q[1],q[2])
yfit2 = func2(xfit2,q[0],q[1])

print('\nBest fit (exponential):')
#print(u'a = {:g} \u00B1 {:g}'.format(q[0],qerr[0]))
#print(u'b = {:g} \u00B1 {:g}'.format(q[1],qerr[1]))
#print(u'c = {:g} \u00B1 {:g}'.format(q[2],qerr[2]))
print(u'b = {:g} \u00B1 {:g}'.format(q[0],qerr[0]))
print(u'c = {:g} \u00B1 {:g}'.format(q[1],qerr[1]))

#-----------------------------------------------------------------------------
# Set up the plot for the fitted curves
#
plt.plot(xfit, yfit, 'k-',linewidth=3.0, label='Fit')
plt.text(0.05 ,145., 'Fit with: $y = a + bx\ e^{-cx}$', fontsize=20)
lam = 1./p[2]*1000.
s = '{0:.2f}'.format(lam)
s = '$\lambda _{q,t} = ' + s + '\ mm $'
plt.text(0.1,130., s, fontsize=20)

plt.plot(xfit2, yfit2, 'b--',linewidth=3.0, label='Fit exp')
plt.text(0.05 ,95., 'Fit with: $y = b\ e^{-cx}$', fontsize=20, color='b')
#lam = 1./q[2]*1000.
lam = 1./q[1]*1000.
s = '{0:.2f}'.format(lam)
s = '$\lambda _{q,t} = ' + s + '\ mm $'
plt.text(0.1,80., s, fontsize=20, color='b')

#-----------------------------------------------------------------------------
# show the final product in one plot
#
fig.subplots_adjust(bottom=0.15, left=0.15)
plt.legend()
plt.show()