# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 16:12:03 2017

@author: bartonjo

Deconvolute the measured Temperature profile from TC data with heat 
flux data from the IRTV.  TC data will be read in from a file and 
fitting will be done by "eye" for now.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

#---------------------------------------------------------------------
# Define parameters
#

# User defined parameters
Feff = 3.15e6 
F = 2.2e6       # const flux on surface [W/m2]
    #***phip = 1.61*F      # transient flux power density [W/m2]
elmtime = 3e-3     # [s] ELM duration
elmf = 50            # "ELM" frequency Hz
finite = 0      # if >0 then finite BCs are to be used
data = 2        # if >0 then use data from file, 1 use TC, 2 plot all
    # data attributes from files
shot = 167408
tcname = 10
pname = 19
tmin = 2.5
tmax = 4.
    # Input the flux intensities min and max to average over (IRTV)
#intermin = 0.5*1e2      # 1 MW/m2 = 100 * 1 W/cm2
intermax = 2.5*1e2
elmmin = 3.1*1e2
#elmmax = 1.3*1e2
    # for plotting
Toffset = 80
toffset = 2.15

# Material constants for pure Tungsten samples
#rho = 19.25e3                   # density [kg/m3]
#cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
#k = 173.                        # thermal conductivity [W/m/K]
# Material constants for Mo
rho = 10.2e3                     # density [kg/m3]
cp = .23e3                       # specific heat [J/kg/K]
k = 138.                         # thermal conductivity [W/m/K]
kappa = k/rho/cp                # thermal diffusivity [m2/s]
d = .75e-2                      # measurement depth [m]
x = d                           # [m] from surface

# numerical parameters and array initialization
num = range(10)                 # summation limit
t = np.linspace(1e-4,3.5,5000)     # seconds
Temp = np.zeros((1, np.size(t)))   # temperature [K]
Telms = np.zeros((1, np.size(t)))

# Calculate other useful parameters
numelms = np.max(t)*elmf        # number of ELMs
for i in t:                     # put the elm duration on t scale
    if (i<=elmtime):
        elmt = i
    else:
        break
tau = [t[i] for i in range(np.size(t)) \
        if np.mod(i+1,np.round(np.size(t)/numelms))==0] #ELM event time
ts = np.zeros((np.size(tau), np.size(t)))

# define useful functions
# calculate the integral of the error function
def ierfc(z):
    return z*erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)

# calculate the shifted time for the elm heat flux
def tshift(t,s,telm):
    shift = [i for i in t if i>=telm]
    z = [np.min(t)/10. for i in t if (i-s)<=0]
    z.extend(shift)
    z = np.array(z[0:np.size(t)])
    return z

# sort of heavyside function
# takes the array ts and returns array of same size with 0s at
# ts indexes equalling mint and 1s otherwise
def heavy(ts,mint):
    h = np.piecewise(ts,[ts==mint, ts>mint],[0,1])
    return h

for m in range(np.size(tau)):
    ts[m][:] = tshift(t,tau[m],elmt)

#---------------------------------------------------------------------
# Calculate ELM hf from Feff and F (inter-ELM)
#

# calculate heavyside summation
Hsum = 0
for m in range(np.size(tau)):       #loop over each elm
    ts2 = heavy(ts[m][:],np.min(t)/10.)*ts[m][:]
    s = 0
    for n in num:
        s += np.exp(-(d*(n+1))**2/kappa/(ts2+1e-10)) + \
                np.exp(-(d*n)**2/kappa/(ts2+1e-10))
    for i in range(np.size(t)):
        if ts2[i]>0:
            if finite>0:
                # finite bcs
                ts2[i]=s[i]/np.sqrt(ts2[i]*np.pi)
            else:
                # semi-infinite bcs
                ts2[i]=1/np.sqrt(ts2[i])
    Hsum += ts2
def moving_av(a, n):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:]-ret[:-n]
    return ret[n-1:]/n    
H = moving_av(Hsum,int(np.size(t)/50.))
end = np.size(H)

# calculate the fraction of the ELM heat flux that it would need to
# be to give an equivalent Temperature profile with a constant hf

if finite>0:
    cc = 0
    for n in num:
        cc += -ierfc((2*d*(n+1))/np.sqrt(4*kappa*t)) - \
                ierfc((2*d*n)/np.sqrt(4*kappa*t))
    dum1 = H/np.sqrt(t[:end])/cc[:end]
    dum2 = np.mean(dum1[int(np.size(dum1)/2.):])
    av_factor = 2*dum2
else:
    dum1 = H/np.sqrt(t[:end]) # this is pretty much a straight line
    dum2 = np.mean(dum1[int(np.size(dum1)/2.):])
    av_factor = elmtime*dum2       # see derivation

# Given the effective hf, Feff, and F, find phip
if finite>0:
    # finite bcs
    phip_av = (Feff-F)/av_factor  # ELM energy density [J/m2]
    # convert phip_av to flux density with finite BCs
    c,s=0,0
    for n in num:
        c += -ierfc(d*(n+1)/np.sqrt(kappa*elmtime)) - \
                ierfc(d*n/np.sqrt(kappa*elmtime))
        s += np.exp(-(d*(n+1))**2/kappa/elmtime) + \
                np.exp(-(d*n)**2/kappa/elmtime)
    phip = phip_av/2/elmtime/np.sqrt(np.pi)*s/c  # flux density W/m2
    phi =  phip_av   
else:
    # semi-infinite bcs
    phip = (Feff-F)/av_factor
    # convert to energy density with semi-infinite BCs
    phi = 2*phip*elmtime

h0 = 'effective hf = {:0.2f} MW/m2'.format(Feff/1e6)
h1 = 'inter-ELM hf = {:0.2f} MW/m2'.format(F/1e6)
h2 = 'ELM hf = {:0.2f} MW/m2'.format(phip/1e6)
print('\nFrom the model:')
print(h0 + '\n' + h1 +'\n' + h2)
print(phip/F)


#---------------------------------------------------------------------
# Evaluate solution of linear heat eaquation with constant and trans
# heat flux on the surface
#

# eval the temp vs time with a constant background heat flux, F
if finite>0:
    # finite BCs
    c = 0
    for n in num:
        c += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*t)) - \
                ierfc((2*d*n+x)/np.sqrt(4*kappa*t))
    Tinter = F/k * np.sqrt(4*kappa*t) * c
else:
    # semi-infinite BCs
    Tinter = F/k * (2*np.sqrt(kappa*t/np.pi)*\
            np.exp(-x**2/4/kappa/t) - x*erfc(x/2/np.sqrt(kappa*t)))

# evaluate the temp vs time for transient ELM heat fluxes, phi
for m in range(np.size(tau)):
    if finite>0:
        # finite BCs
        s=0
        for n in num:
            s += np.exp(-(2*d*(n+1)-x)**2/4/kappa/ts[m][:]) + \
                    np.exp(-(2*d*n+x)**2/4/kappa/ts[m][:])
        Telms += phi/k*np.sqrt(kappa/np.pi/ts[m][:])* s * \
                np.piecewise(ts[m][:],
                [ts[m][:]==np.min(t)/10., ts[m][:]>np.min(t)/10.],
                [0,1])
    else:
        # semi-infinite BCs
        Telms += phi/k*np.sqrt(kappa/np.pi/ts[m][:]) * \
                np.exp(-x**2/4/kappa/ts[m][:]) * \
                heavy(ts[m][:],np.min(t)/10.)

# Linearity allows us to add the inter-ELM and ELM solutions
Temp = Tinter + Telms

# evaluate the effective temperature with Feff
if finite>0:
    # finite BCs
    c = 0
    for n in num:
        c += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*t)) - \
                ierfc((2*d*n+x)/np.sqrt(4*kappa*t))
    Teff = Feff/k * np.sqrt(4*kappa*t) * c
else:
    # semi-infinite BCs
    Teff = Feff/k * (2*np.sqrt(kappa*t/np.pi)*np.exp(-x**2/4/kappa/t) - \
            x*erfc(x/2/np.sqrt(kappa*t)))

#---------------------------------------------------------------------
# Import TC data and anything else from files
#
if data>0:
    # assemble file names
    ftc = 'decondata/'+'hfdecon_'+str(shot)+'_'+str(tcname)+'tc.dat'
    flp = 'decondata/'+'hfdecon_'+str(shot)+'_'+str(pname)+'p.dat'
    fir = 'decondata/'+'hfdecon_'+str(shot)+'_ir.dat'
    fe = 'decondata/'+'hfdecon_'+str(shot)+'_elmf.dat'
    fts = 'decondata/'+'hfdecon_'+str(shot)+'_dts.dat'
    ffs = 'decondata/'+'hfdecon_'+str(shot)+'_fs5.dat'
    flpmed = 'decondata/'+'hfdecon_'+str(shot)+'_'+\
            str(pname)+'pmed.dat'    
    
    # get data from file and put into np arrays
    tct, tcTemp, tcderiv = np.loadtxt(ftc, unpack=True)
    lpt, jsat, Te, ang, lpq = np.loadtxt(flp, unpack=True)
    irt, irq = np.loadtxt(fir, unpack=True)
    dsselmt, dsselmf = np.loadtxt(fe, unpack=True)
    tst, tsq, tsq1 = np.loadtxt(fts, unpack=True)
    fst, fs = np.loadtxt(ffs, unpack=True)
    lptmed, lpqmed = np.loadtxt(flpmed, unpack=True)
    
    # for ease of plot manipulation, decrease the array size
    tct1 = [tct[i] for i in range(np.size(tct)) if (i+1)%100==0]
    tcTemp1 = [tcTemp[i] for i in range(np.size(tct)) if (i+1)%100==0]
    tct1 = np.array(tct1)
    tcTemp1 = np.array(tcTemp1)

    # calculate ELM and inter-ELM heat fluxes from IRTV and print
        # Find time region of interest
    zmin = np.min(np.where(irt/1e3>tmin))
    zmax = np.min(np.where(irt/1e3>tmax))
        # Sometimes may need to subtract out background slope
    def sub_line(x1,y1,x2,y2):
        m = (y2-y1)/(x2-x1)
        b = y2-m*x2
        return m,b
    m, b = sub_line(irt[zmin]/1e3,irq[zmin]/1e2,
                    irt[zmax]/1e3,irq[zmax]/1e2)
    y = m*irt/1e3+b
    back = 1.6      # "flat" background to add back in
        # find the averages and print them [MW/m2]
    interarr = [irq[i] for i in range(np.size(irq[zmin:zmax]))+zmin\
            if (irq[i]<intermax)]  
#    interarr = [irq[i]-y[i]*1e2+back*1e2\
#            for i in range(np.size(irq[zmin:zmax]))+zmin\
#            if (irq[i]-y[i]*1e2+back*1e2<intermax)]
    interarr = np.array(interarr)
    elmarr = [irq[i] for i in range(np.size(irq[zmin:zmax]))+zmin\
            if (irq[i]>elmmin)]
#    elmarr = [irq[i]-y[i]*1e2+back*1e2\
#            for i in range(np.size(irq[zmin:zmax]))+zmin\
#            if (irq[i]-y[i]*1e2+back*1e2>elmmin)]
    elmarr = np.array(elmarr)
    interave = np.mean(interarr)/1e2    # MW/m2
    elmave = np.mean(elmarr)/1e2
    hinter = '{:0.2f}'.format(interave)
    helm = '{:0.2f}'.format(elmave)
    print('\nFrom the data:')
    print('IRTV inter-ELM hf = '+hinter+' MW/m2')
    print('IRTV ELM hf = '+helm+' MW/m2')
    print(elmave/interave)
        # plot results
    if data==2:
        plt.subplots()
        plt.plot(irt/1e3,irq/1e2)
        plt.plot(irt[zmin:zmax]/1e3,irq[zmin:zmax]/1e2)
#        plt.plot([0,6],[intermin/1e2,intermin/1e2],'c', linewidth=3)
        plt.plot([0,6],[intermax/1e2,intermax/1e2],'c', linewidth=3)
        plt.plot([0,6],[elmmin/1e2,elmmin/1e2],'m', linewidth=3)
#        plt.plot([0,6],[elmmax/1e2,elmmax/1e2],'m', linewidth=3)
        plt.plot([0,6],[interave,interave],'r')
        plt.plot([0,6],[elmave,elmave],'r')
#        plt.plot(irt/1e3,y)
#        plt.plot(irt[zmin:zmax]/1e3,
#                 irq[zmin:zmax]/1e2-y[zmin:zmax]+back)
        plt.title('IRTV')
        plt.xlabel('t [s]', fontsize=18)
        plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=18)
        plt.show()
    
    # convert DTS q_parallel to q_perp
    zmin = np.min(np.where(lpt/1e3>tmin))
    zmax = np.min(np.where(lpt/1e3>tmax))
    av_ang = np.mean(ang[zmin:zmax])
    tsq = tsq*np.sin(np.radians(av_ang))
    tsq1 = tsq1*np.sin(np.radians(av_ang))
    
    # Calculate inter-ELM hf from LPs
#    lpqmax = 30.0e2      # W/cm2
#    lpqarr = [lpq[i] for i in range(np.size(lpq[zmin:zmax]))+zmin\
#            if (lpq[i]<lpqmax)]
#    lpqarr = np.array(lpqarr)
    lpqarr = lpq[zmin:zmax]
    lpqave = np.mean(lpqarr)/1e2
    hlpqave = '{:0.2f}'.format(lpqave)
    print('\nLP inter-ELM hf = '+hlpqave+' MW/m2')
    if data==2:
        plt.subplots()
        plt.plot(lpt/1e3,lpq/1e2)
        plt.plot(lpt[zmin:zmax]/1e3,lpq[zmin:zmax]/1e2)
        plt.plot([0,6],[lpqave,lpqave],linewidth=2,label='raw ave')
        
        zmin = np.min(np.where(lptmed/1e3>tmin))
        zmax = np.min(np.where(lptmed/1e3>tmax))
    #    lpqarr=[lpqmed[i] for i in range(np.size(lpqmed[zmin:zmax]))+zmin\
    #            if (lpqmed[i]<lpqmax)]
    #    lpqarr = np.array(lpqarr)
        lpqarr = lpqmed[zmin:zmax]
        lpqave = np.mean(lpqarr)/1e2
        hlpqave = '{:0.2f}'.format(lpqave)
        print('LP inter-ELM hf = '+hlpqave+' MW/m2 (median filter)')
        plt.plot(lptmed/1e3,lpqmed/1e2,'s',markersize=10)
        plt.plot(lptmed[zmin:zmax]/1e3,lpqmed[zmin:zmax]/1e2,'s',
                 markersize=10)
        plt.plot([0,6],[lpqave,lpqave],linewidth=2,label='filtered ave')
#        plt.ylim([0,10])
        plt.title('LP heat flux')
        plt.xlabel('t [s]', fontsize=18)
        plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=18)
        plt.legend(loc='best')
        plt.show()

#---------------------------------------------------------------------
# Plotting
#

# Calculated Temperature vs time
fig, ax = plt.subplots()

h = '{:0.2f}'.format(d*100)
plt.plot(t+toffset, Teff+Toffset, linewidth=3, 
         label=r'$T_{eff}$ (' +h+ ' cm depth)')
plt.plot(t+toffset, Temp[0][:]+Toffset, '--', linewidth=3, 
         label=r'$T_{inter}+T_{ELMs}$ ('+h+' cm depth)')
plt.plot(t+toffset, Tinter+Toffset, linewidth=3, 
         label=r'$T_{inter}$')
plt.plot(t+toffset, Telms[0][:]+Toffset, linewidth=3, 
         label=r'$T_{ELMs}$')
if data>0:
    plt.plot(tct1/1000., tcTemp1, 'c', label='data')
    
# To illustrate <Telms>:
#dum = H*phi/k*np.sqrt(kappa/np.pi)
#dumt = t[:end]
#x = dumt+toffset
#y = dum+Toffset
#plt.plot(x[::10],y[::10],'r', linewidth = 3, label=r'$<T_{ELMs}>$')
#plt.title('@ x = 0', fontsize=16)
#plt.xlim([toffset, 5.5])
#plt.ylim([Toffset, 200])

plt.plot([0, np.max(t)+2], [0,0], 'k')
plt.xlim([1, 6])
plt.ylim([0, 500])

h1 = r'$F_{inter}=$' 
h = ' {:0.1f} MW/m$^2$ and\n'.format(F/1e6)
h1 = h1 + h
h2 = r'$\phi_{ELM}=$'
h = ' {:0.1f}'.format(phip/1e6)
h2 = h2 + h
h3 = ' MW/m$^2$ at {:0.0f} Hz'.format(elmf)
h4 = r'$F_{eff}=$'
h = ' {:0.1f} MW/m$^2$'.format(Feff/1e6)
plt.title(h1+h2+h3, fontsize=16)
#plt.title(h4+h, fontsize=16)
plt.xlabel('Time [s]', fontsize=18)
plt.tick_params(axis='x', labelsize=16)
plt.ylabel('Temperature [C]', fontsize=18)
plt.tick_params(axis='y', labelsize=16)
plt.legend(loc='upper left')
fig.subplots_adjust(top=0.89)
plt.show()















