# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 16:12:03 2017

@author: bartonjo

Deconvolute the inter-ELM and ELM contributions to the measured 
Temperature profile from TC data. Heat fluxes from TC and LP 
calculations will be compared to heat flux data from the IRTV.  All 
data will be read in from a file.  There is only one free parameter: 
the peak ELM heat flux.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import import_hfdecon_data as decondata
from scipy.optimize import curve_fit
#from scipy.stats import chisquare
from scipy.stats import pearsonr
import h5py

#---------------------------------------------------------------------
# User defined parameters
#

# hdf5 file name
filename = 'filename.hdf5'

# data attributes from files
shot = 167297
tcname = 9      # 9, 10, 11 for metal rings
pname = 23      # 51, 19, (21), 15, (17)
tmin = 2.07     
tmax = 4.5      # tmin and tmax is the time window for averaging data 

# choose BCs
finite = 1      # 0 semi-infinite, 1 finite BCs

# Thickness of slab, d, and location of T analysis, x
d = 6e-3#.75e-2                      # measurement depth [m]
angle = 0                       # angle of button from horizontal
# .75e-2 insert TC, 2.77e-3 (3.18) 15o angled button, 1.57e-3 button
if angle>0:
    d += .75*6e-3*np.tan(np.radians(angle))
x = 0                           # [m] from surface

# force the values of the heat flux and do no fitting if =1
force = 0
#Fmax = 1.24e6 #infinite 23.84e6  #finite 1.24
#Feff = 1.56e6 #infinite 3.16e6 #finite 1.56

# Material constants for pure Tungsten samples
rho = 19.25e3                   # density [kg/m3]
cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
k = 125.                        # thermal conductivity [W/m/K]
                                # 173 RT, 125 (500C), 112 (927C)
# Material constants for Mo
#rho = 10.2e3                     # density [kg/m3]
#cp = .23e3                       # specific heat [J/kg/K]
#k = 138.                         # thermal conductivity [W/m/K]
# thermal diffusivity [m2/s]
kappa = k/rho/cp                

# other numerical parameters
if force==0:
    # times with constant Finter, =0 means const always
    interstep = 0.  
else:
    interstep = 0.
num = range(20)                 # summation limit for finite BCs
t = np.linspace(1e-4,3.5,5000)     # seconds

#---------------------------------------------------------------------
# Extract useful parameters from the data
#

s = decondata.getdecondata(shot=shot, tcname=tcname, pname=pname)

# return array with ELMs separated out and array with only ELMs
def elm_sep(t,dat,s=s):
    """(t, dat) dat is data array with corresponding time array, t"""
        # separate elms in ir data
    for n in range(np.size(s.elmstart)):
        t1 = s.elmstart[n]  # [ms]
        t2 = t1+s.dur[n]
        # find elements in t when this ELM occurs
        try:
            zelm1 = np.min(np.where(t>=t1))
            zelm2 = np.min(np.where(t>=t2))
        except ValueError:
#            print(np.max(t))
#            print(t1)
            break
        if n==0:
            zelms = list(np.arange(zelm1,zelm2+1))
        else:
            zelms += list(np.arange(zelm1,zelm2+1))
    elmt = [t[i] for i in zelms]
    elmd = [dat[i] for i in zelms]
    intert = [t[i] for i in range(np.size(t)) if i not in zelms]
    interd = [dat[i] for i in range(np.size(t)) if i not in zelms]
    return intert, interd, elmt, elmd

# parse in time the inter-ELM heat flux data and find the hf in those
# parsed time windows. Return the time window data and the hfs.
def parse_inter(t,dat,interstep=interstep):
    '''(t,dat) assummed to be within tmin and tmax. t in [s]!!!'''
    t = t[:np.min(np.where(t==np.max(t)))]
    dat = dat[:np.min(np.where(t==np.max(t)))]
    # number of slices
    num = np.int(np.ceil((np.max(t)-np.min(t))/interstep))
    # find where these slices begin
    steps = [t[0]+interstep*i for i in range(num)]
    steps = np.array(steps)
    tstart = [t[np.min(np.where(t>steps[i]))] for i in range(num)]    
    tstart = np.array(tstart)
    # slices end one element behind where they begin
    tend=[t[np.min(np.where(t==tstart[i+1]))-1] for i in range(num-1)]
    tend += [t[-1]]
    tend = np.array(tend)
    # find the heat flux within the slices
    F = np.zeros(np.size(tstart))
    for i in range(num):
        zmin = np.min(np.where(t==tstart[i]))
        zmax = np.min(np.where(t==tend[i]))
#        print(zmin)
#        print(zmax)
        F[i] = np.mean(dat[zmin:zmax+1])
#    print(F)
    return tstart, tend, F

# Find the time offset when the strike point is on the TC
tcloc = {9:1.41,10:1.354,11:1.324}
osp = np.array([s.osp[i] for i in range(np.size(s.osp)) \
            if s.osp[i]>0 and s.ospt[i]>2000.])
ospt = [s.ospt[i] for i in range(np.size(s.osp)) \
            if s.osp[i]>0 and s.ospt[i]>2000.]
ospt = np.array(ospt)
z1 = np.min(np.where(osp<tcloc[tcname+0]))
z2 = np.min(np.where(osp[(z1+1):]<tcloc[tcname+0]))
z3 = 1+z1+z2
toffset = ospt[z3]/1e3  # convert ms to s
toffset = toffset + t[0]   # simulation starts at t[0] not zero

# Find the wetted area of the ring
# average OSP radius
z4 = np.max(np.where(ospt/1e3<tmax))
r = 1.41#np.mean(osp[z3:z4+1])
if r>1.4:  # on the shelf ring
    A = np.pi * (1.45**2 - r**2)
    print('Getting area with OSP on shelf ring')
    print('{:0.2f}'.format(A)+' m2')
elif (r<1.4 and r>1.37) or tcname == 9:  # shelf ring
    A = np.pi * (1.45**2 - 1.4**2)
    print('Getting area of entire shelf ring')
    print('{:0.2f}'.format(A)+' m2')
elif r>1.32:  # on the floor ring
    A = np.pi * (1.37**2 - r**2)
    print('Getting area with OSP on floor ring')
    print('{:0.2f}'.format(A)+' m2')
else:  # floor ring
    A = np.pi * (1.37**2 - 1.32**2)
    print('Getting area of entire floor ring')
    print('{:0.2f}'.format(A)+' m2')

# for ease of plot manipulation and fitting, decrease the TC array size
tct1 = s.tct[0::100]
tcTemp1 = s.tcTemp[0::100]

# Find the temperature offset (i.e. starting temperature)
z2 = np.min(np.where(tct1/1e3>toffset))
z1 = z2-20
Toffset = np.mean(tcTemp1[z1:z2+1])

# get the inter-ELM heat flux data
zmin = np.min(np.where(s.lpt/1e3>tmin))
zmax = np.min(np.where(s.lpt/1e3>tmax))
av_ang = np.mean(s.ang[zmin:zmax+1])
    # convert DTS q_parallel to q_perp [MW/m2]
tsq = s.tsq*np.sin(np.radians(av_ang))/1e2
tsq1 = s.tsq1*np.sin(np.radians(av_ang))/1e2
    # raw average of LP data
lpqarr = s.lpq[zmin:zmax+1]
lpqave = np.mean(lpqarr)/1e2    # converting W/cm2 to MW/m2
hlpqave = '{:0.2f}'.format(lpqave)
print('\nLP mean hf = '+hlpqave+' MW/m2 (raw average)')
lpt,lpq,let,leq = elm_sep(s.lpt[zmin:zmax+1],s.lpq[zmin:zmax+1])
lpt = np.array(lpt)
lpq = np.array(lpq)
let = np.array(let)
leq = np.array(leq)
lptf,lpqf,letf,leqf = elm_sep(s.lpt,s.lpq)
lptf = np.array(lptf)
lpqf = np.array(lpqf)
letf = np.array(letf)
leqf = np.array(leqf)
# heat flux used in calculating inter-ELM temperature history
lpqinter = np.mean(lpq)/1e2
hlpqave = '{:0.2f}'.format(lpqinter)
print('LP inter-ELM hf = '+hlpqave+' MW/m2 (ELM data removed)')
    # median-filter average of LP data
zmin = np.min(np.where(s.lptmed/1e3>tmin))
zmax = np.min(np.where(s.lptmed/1e3>tmax))
lpqarr = s.lpqmed[zmin:zmax+1]
lpqavemed = np.mean(lpqarr)/1e2
hlpqave = '{:0.2f}'.format(lpqavemed)
print('LP filtered hf = '+hlpqave+' MW/m2 (median filter)')
if interstep > 0:
    t1, t2, Finter = parse_inter(lptf/1e3,lpqf/1e2*1e6)
else:
    t1,t2=0,0
    Finter = lpqinter*1e6     # converted to W/m2

# get ELM properties
zmin = np.min(np.where(s.elmstart/1e3>toffset))
    #simulation runs for t[-1] seconds from when OSP is on target
#try:
#    zmax = np.min(np.where(s.elmstart/1e3>(tmax)))
#except ValueError:
#    zmax = np.min(np.where(s.elmstart/1e3>tmax))
zmax = np.size(s.elmstart) - 1
    # ELM events
tau = s.elmstart[zmin:zmax+1]/1e3   # seconds
    # ELM duration
dtau = s.dur[zmin:zmax+1]/1e3       # seconds
    # put taustart and tauend on simulation (i.e. t) scale
taustart = tau - toffset
tauend = (tau+dtau) - toffset
    # find relative energy fraction of each ELM
emax = np.max(s.esize[zmin:zmax+1])
efrac = s.esize[zmin:zmax+1]/emax

#---------------------------------------------------------------------
# Define functions for model
#

# calculate the integral of the error function
def ierfc(z):
    return z*erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)
    
# calculate the heavyside function
def hside(c):
    return np.piecewise(c, [c>=0, c<0], [1,0])

# calculate the temperature due to a single constant heat flux
def T_hf(F,x=x,t=t,k=k,kappa=kappa,finite=finite,num=num,d=d):
    if finite>0:
        c = 0
        for n in num:
            c += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*t)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*t))      
        T = F/k * np.sqrt(4*kappa*t) * c
    else:
        T = F/k * np.sqrt(4*kappa*t) * (-ierfc(x/np.sqrt(4*kappa*t)))
    return T

# return t minus ELM start and stop time (dt1 and dt2)
def tshift(tstart,tend,t=t):
    dt1 = t - tstart
    dt2 = t - tend
    #remove zeros from dt, if any, to avoid dividing by 0
    y1 = [i for i in range(np.size(dt1)) if dt1[i]==0.]
    y2 = [i for i in range(np.size(dt2)) if dt2[i]==0.]
    dt1 = np.array([dt1[i] if i not in y1 else 0\
            for i in range(np.size(dt1))])
    dt2 = np.array([dt2[i] if i not in y2 else 0\
            for i in range(np.size(dt2))])
    return dt1, dt2

# calculate the temperature due to both inter-ELM and ELM heat flux
def T_tot(Fmax,Finter=Finter,t=t,t1=t1,t2=t2,toffset=toffset,
          efrac=efrac,taustart=taustart,tauend=tauend):
    #1-d heat equation temperature history inter-ELM term
    if np.size(Finter) == 1:  # constant inter-ELM hf
        Tinter = T_hf(Finter,t=t)
    else:
        Tinter = np.zeros((np.size(t)))
        # accomodate shorter t ranges
#        zmin = np.min(np.where((t1-toffset)>np.min(t)))
        zmax = np.max(np.where((t1-toffset)<np.max(t)))
        for i in np.arange(0,zmax+1):
            dt1, dt2 = tshift((t1[i]-toffset),(t2[i]-toffset),t=t)
            #1-d heat equation temperature history ELM terms
            Tsum1 = T_hf(Finter[i],t=np.abs(dt1))*hside(dt1)
            Tsum2 = T_hf(Finter[i],t=np.abs(dt2))*hside(dt2)
            Tinter += (Tsum1 - Tsum2)
    #temperature due to all ELMs
    Telms = np.zeros((np.size(Tinter)))
    # accomodate shorter t ranges
#    zmin = np.min(np.where(taustart>np.min(t)))
    zmax = np.max(np.where(taustart<np.max(t)))
    for m in np.arange(0,zmax+1):
        dt1, dt2 = tshift(taustart[m],tauend[m],t=t)
        #1-d heat equation temperature history ELM terms
        Tsum1 = T_hf(Fmax*efrac[m],t=np.abs(dt1))*hside(dt1)
        Tsum2 = T_hf(Fmax*efrac[m],t=np.abs(dt2))*hside(dt2)
        Telms += (Tsum1 - Tsum2)
    T = Tinter + Telms
    return T

#---------------------------------------------------------------------
# Fit TC time history using a single constant (effective) heat flux
# 

if force == 0:
    # define the region of the data that you want to fit (tmin,tmax)
    zmin = np.min(np.where(tct1/1e3>tmin))
    zmax = np.min(np.where(tct1/1e3>tmax))
    tctdata = tct1[zmin:zmax+1]/1e3-toffset
    tcdata = tcTemp1[zmin:zmax+1]-Toffset
    
    # define the function in the way curve_fit likes: f(xdata, *params)
    def single_hf_fit(time, F):
        return T_hf(F,t=time)
    
    # fit the data
    popt, pcov = curve_fit(
            single_hf_fit, tctdata, tcdata, bounds=(1e5,1e9))
    Feff = popt[0]  # change from array to float    
    
    # get the std dev error of the fit
    perr = np.max(np.sqrt(np.diag(pcov)))
    
    # get the chi square of the fit
    fit = T_hf(Feff,t=tctdata)
#    chiv, p = chisquare(tcdata,f_exp = fit, ddof=1)
#    print(chiv) 
    rsq = pearsonr(tcdata,fit)[0]
    
    # calculate the temperature history for all t
    Teff = T_hf(Feff)
    
    # calculate the deposited energy
    effflux = Feff * (tmax-tmin)    # J/m2
    
    # print out results
    hperr = '{:0.4f}'.format(perr/1e6)
    hFeff = '{:0.2f}'.format(Feff/1e6)
    hefff = '{:0.2f}'.format(effflux/1e6)
    print('\nFeff = '+hFeff+' MW/m2')
    print('std_dev of Feff = '+hperr+' MW/m2')
#    p = '{:0.6f}'.format(p)    
#    print('Chi-square test of fit = '+p)
    r = '{:0.6f}'.format(rsq)
    print('linear R-sq of fit = '+r)
    print('Effective energy = '+hefff+' MJ/m2')
else:
    # find equivalent perp heat flux for a perscribed angle
    if angle>0:
        a = angle+av_ang
        factor = np.sin(np.radians(a))/np.sin(np.radians(av_ang))
        Feff = Feff * factor   
    Teff = T_hf(Feff)

#---------------------------------------------------------------------
# Fit TC time history using the inter-ELM heat flux from the data and
# using the ELM heat flux as the free parameter
#

if force==0:
    # tctdata and tcdata is the same as in the above section
    
    # define the function in the way curve_fit likes: f(xdata, *params)
    def ELM_hf_fit(time, F):
        return T_tot(F,t=time)
    
    # fit the data
    try:
        popt, pcov = curve_fit(ELM_hf_fit, 
                           tctdata, tcdata, bounds=(1e4,1e9))
    except ValueError:
        popt, pcov = curve_fit(ELM_hf_fit, 
                           tctdata, tcdata)
    Fmax = popt[0]  # change from array to float
    
    # get the std dev error of the fit
    perr = np.max(np.sqrt(np.diag(pcov)))
    
    # get the chi square of the fit
#    chiv, p = chisquare(tcdata,f_exp = T_tot(Fmax,t=tctdata))
#    print(chiv) 
    rsq = pearsonr(tcdata,T_tot(Fmax,t=tctdata))[0]
    
    # calculate the temperature history for all t
    # This is Fmax in the fitting time range. Find the real Fmax
#    zmin = np.min(np.where(taustart>np.min(tctdata)))
#    zmax = np.max(np.where(taustart<np.max(tctdata)))
#    Fmax = Fmax/np.max(efrac[zmin:zmax])
    Tfit = T_tot(Fmax)
    
    # calculate the deposited energy
    z = np.max(np.where(tau<tmax))
    elmflux = [Fmax*efrac[m]*dtau[m] for m in range(np.size(dtau[:z]))]
    elmflux = np.array(elmflux)         # J/m2  
    elmenflux = np.sum(elmflux)         # J/m2
    interflux = np.mean(Finter) * (tmax-tmin)    # J/m2
    
    # print out results
    hperr = '{:0.4f}'.format(perr/1e6)
    hFmax = '{:0.2f}'.format(Fmax/1e6)
    hFave = '{:0.2f}'.format(Fmax/1e6*np.mean(efrac))
    helmf = '{:0.2f}'.format(elmenflux/1e6)
    hintf = '{:0.2f}'.format(interflux/1e6)
    htot = '{:0.2f}'.format((elmenflux+interflux)/1e6)
    print('\nFmax = '+hFmax+' MW/m2')
    print('std_dev of Fmax = '+hperr+' MW/m2')
#    p = '{:0.6f}'.format(p)    
#    print('Chi-square test of fit = '+p)
    r = '{:0.6f}'.format(rsq)
    print('linear R-sq of fit = '+r)
    print('Felm_ave = '+hFave+' MW/m2')
    print('ELM energy = '+helmf+' MJ/m2')
    print('inter-ELM energy = '+hintf+' MJ/m2')
    print('total energy = '+htot+' MJ/m2')
else:
    # find equivalent perp heat flux for a perscribed angle
    if angle>0:
        a = angle+av_ang
        factor = np.sin(np.radians(a))/np.sin(np.radians(av_ang))
        Finter = Finter * factor
        Fmax = Fmax * factor
    Tfit = T_tot(Fmax,Finter=Finter)
    # calculate the deposited energy
    z = np.max(np.where(tau<tmax))
    elmflux = [Fmax*efrac[m]*dtau[m] for m in range(np.size(dtau[:z]))]
    elmflux = np.array(elmflux)         # J/m2  
    elmenflux = np.sum(elmflux)         # J/m2
    interflux = Finter * (tmax-tmin)    # J/m2

#---------------------------------------------------------------------
# Plotting
#

# Calculated Temperature vs time
fig, (ax1,ax2) = plt.subplots(2,sharex=True)

#h = '{:0.2f}'.format(x*100)
ax1.plot(t+toffset, Tfit+Toffset, linewidth=3, 
         label=r'$T(F_{inter},F_{ELMs})$')# ('+h+' cm depth)')
ax1.plot(t+toffset, Teff+Toffset, '--', linewidth=1,
         label=r'$T(F_{eff})$')# (' +h+ ' cm depth)')
ax1.plot(tct1/1000., tcTemp1, 'k', label='data')

F = np.zeros((np.size(t)))
if np.size(Finter)>1:    
    j=0
    for i in range(np.size(F)):
        if (t[i]+toffset)>t2[j]:
            j+=1
        if j>=np.size(Finter):
            F[i:] = Finter[j-1]
            break
        else:
            F[i]=Finter[j]
else:
    F = F + Finter
for m in range(np.size(taustart)):
    felm = (hside(t-taustart[m])-hside(t-tauend[m]))*efrac[m]*Fmax
    F += felm
ax2.plot(t+toffset,F/1e6,
         label=r'$F_{inter} + F_{ELMs}$')
ax2.plot([t[0]+toffset,t[-1]+toffset],[Feff/1e6,Feff/1e6],'--',
         label=r'$F_{eff}$')

#plt.plot([0, np.max(t)+2], [0,0], 'k')
#plt.xlim([0, 3.5])
#plt.ylim([0, 300])

h = '{:0.2f}'.format(x*100)
sh = '{:0.0f}'.format(shot)
ax1.set_title('shot '+sh+'\nTemperature at '+h+' cm depth',
              fontsize=16)
#h1 = r'$F_{inter}=$' 
#h = ' {:0.2f} MW/m$^2$ and '.format(np.mean(Finter)/1e6)
#h1 = h1 + h
#h2 = r'$<F_{ELM}>=$'
#h = ' {:0.2f}'.format(Fmax/1e6*np.mean(efrac))
#h2 = h2 + h + ' MW/m$^2$'
#h3 = r'$F_{eff}=$'
#h = ' {:0.2f} MW/m$^2$'.format(Feff/1e6)
#h3 = h3 + h
#plt.title(h1+h2+'\n'+h3, fontsize=16)
#ax1.set_title(h1+h2, fontsize=16)
#plt.title(h4+h, fontsize=16)
ax2.set_xlabel('Time [s]', fontsize=18)
ax2.tick_params(axis='x', labelsize=16)
ax1.set_ylabel('Temperature [C]', fontsize=18)
ax2.set_ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=18)
ax1.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax1.legend(loc='upper left')
ax2.legend(loc='best')
#fig.subplots_adjust(top=0.89)
fig.subplots_adjust(hspace=0.1)
plt.show()



#---------------------------------------------------------------------
# Define plotting functions that can be called after running
# this script
# 

# Temperature due to just ELMs
def ELMtemp(t=(t+toffset),numELMs=np.size(taustart),taustart=taustart,
            tauend=tauend,efrac=efrac,Fmax=Fmax):
    Telms = np.zeros((np.size(t)))
    for m in range(numELMs):
        dt1, dt2 = tshift(taustart[m],tauend[m])
        #1-d heat equation temperature history ELM terms
        Tsum1 = T_hf(Fmax*efrac[m],t=np.abs(dt1))*hside(dt1)
        Tsum2 = T_hf(Fmax*efrac[m],t=np.abs(dt2))*hside(dt2)
        Telms += (Tsum1 - Tsum2)
    plt.subplots()
    plt.plot(t,Telms)
    plt.title('Telms')

# heat flux history
def hfhist(t=t,Finter=Finter,taustart=taustart,tauend=tauend,
           efrac=efrac,Fmax=Fmax,toffset=toffset):
    '''Assumes Finter is scalar!'''
    F = Finter+np.zeros((np.size(t)))
    for m in range(np.size(taustart)):
        felm = (hside(t-taustart[m])-hside(t-tauend[m]))*efrac[m]*Fmax
        F += felm
    plt.subplots()
    plt.plot(t+toffset,F/1e6)
    return t+toffset,F/1e6

# return coefficients for line, given two points
def sub_line(x1,y1,x2,y2):
        m = (y2-y1)/(x2-x1)
        b = y2-m*x2
        return m,b

# IR heat flux plotting !!!!! should be integrating over the width of 
# the heat flux wetted area of the insert!!!!!!
def irplot(s=s,tmin=tmin,tmax=tmax):
    # calculate ELM and inter-ELM heat fluxes from IRTV and print
        # Find time region of interest
    zmin = np.min(np.where(s.irt/1e3>tmin))
    zmax = np.min(np.where(s.irt/1e3>tmax))
        # Sometimes may need to subtract out background slope
#    m, b = sub_line(s.irt[zmin]/1e3,s.irq[zmin]/1e2,
#                    s.irt[zmax]/1e3,s.irq[zmax]/1e2)
#    y = m*s.irt/1e3+b
#    back = 1.6      # "flat" background to add back in
        # find the averages and print them [MW/m2]
#    interarr = [irq[i] for i in range(np.size(irq[zmin:zmax]))+zmin\
#            if (irq[i]<intermax)]  
##    interarr = [irq[i]-y[i]*1e2+back*1e2\
##            for i in range(np.size(irq[zmin:zmax]))+zmin\
##            if (irq[i]-y[i]*1e2+back*1e2<intermax)]
#    interarr = np.array(interarr)
#    elmarr = [irq[i] for i in range(np.size(irq[zmin:zmax]))+zmin\
#            if (irq[i]>elmmin)]
##    elmarr = [irq[i]-y[i]*1e2+back*1e2\
##            for i in range(np.size(irq[zmin:zmax]))+zmin\
##            if (irq[i]-y[i]*1e2+back*1e2>elmmin)]
#    elmarr = np.array(elmarr)
#    interave = np.mean(interarr)/1e2    # MW/m2
#    elmave = np.mean(elmarr)/1e2
#    hinter = '{:0.2f}'.format(interave)
#    helm = '{:0.2f}'.format(elmave)
#    print('\nFrom the data:')
#    print('IRTV inter-ELM hf = '+hinter+' MW/m2')
#    print('IRTV ELM hf = '+helm+' MW/m2')
#    print(elmave/interave)
        # plot results

    plt.subplots()
    plt.plot(s.irt/1e3,s.irq/1e2)
    plt.plot(s.irt[zmin:zmax+1]/1e3,s.irq[zmin:zmax+1]/1e2)
    # print out the energy deposited in this time window
    dt = s.irt[zmin+10]/1e3 - s.irt[zmin+9]/1e3
    e = np.sum(s.irq[zmin:zmax+1])*1e4*dt
    henergy = '{:0.2f}'.format(e/1e6)
    print('\nIR deposited energy = '+henergy+' MJ/m2')
##        plt.plot([0,6],[intermin/1e2,intermin/1e2],'c', linewidth=3)
#        plt.plot([0,6],[intermax/1e2,intermax/1e2],'c', linewidth=3)
#        plt.plot([0,6],[elmmin/1e2,elmmin/1e2],'m', linewidth=3)
##        plt.plot([0,6],[elmmax/1e2,elmmax/1e2],'m', linewidth=3)
#        plt.plot([0,6],[interave,interave],'r')
#        plt.plot([0,6],[elmave,elmave],'r')
##        plt.plot(irt/1e3,y)
##        plt.plot(irt[zmin:zmax]/1e3,
##                 irq[zmin:zmax]/1e2-y[zmin:zmax]+back)
    plt.title('IRTV')
    plt.xlabel('t [s]', fontsize=18)
    plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=18)
    plt.show()
    irt,irq,iet,ieq = elm_sep(s.irt[zmin:zmax+1], s.irq[zmin:zmax+1])
    irqave = np.mean(np.array(irq))/1e2    # converting W/cm2 to MW/m2
    hirqave = '{:0.2f}'.format(irqave)
    print('\nIR inter-ELM hf = '+hirqave+' MW/m2')
    ieqave = np.mean(np.array(ieq))/1e2    # converting W/cm2 to MW/m2
    hieqave = '{:0.2f}'.format(ieqave-irqave)
    print('IR ELM hf = '+hieqave+' MW/m2')
    return np.array(irt), np.array(irq), np.array(iet), np.array(ieq)
    
    # Calculate inter-ELM hf from LPs
#    lpqmax = 30.0e2      # W/cm2
#    lpqarr = [lpq[i] for i in range(np.size(lpq[zmin:zmax]))+zmin\
#            if (lpq[i]<lpqmax)]
#    lpqarr = np.array(lpqarr)

def lpplot(s=s,lpqave=lpqave,lpqavemed=lpqavemed,tmin=tmin,tmax=tmax):
    zmin = np.min(np.where(s.lpt/1e3>tmin))
    zmax = np.min(np.where(s.lpt/1e3>tmax))
    plt.subplots()
    plt.plot(s.lpt/1e3,s.lpq/1e2)
    plt.plot(s.lpt[zmin:zmax]/1e3,s.lpq[zmin:zmax]/1e2)
    plt.plot([0,6],[lpqave,lpqave],linewidth=2,label='raw ave')
    
   
#    lpqarr=[lpqmed[i] for i in range(np.size(lpqmed[zmin:zmax]))+zmin\
#            if (lpqmed[i]<lpqmax)]
#    lpqarr = np.array(lpqarr)
    
    plt.plot(s.lptmed/1e3,s.lpqmed/1e2,'s',markersize=10)
    plt.plot(s.lptmed[zmin:zmax]/1e3,s.lpqmed[zmin:zmax]/1e2,'s',
             markersize=10)
    plt.plot([0,6],[lpqavemed,lpqavemed],linewidth=2,
             label='filtered ave')
#        plt.ylim([0,10])
    plt.title('LP heat flux')
    plt.xlabel('t [s]', fontsize=18)
    plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=18)
    plt.legend(loc='best')
    plt.show()

def showlpparse(s=s,tmin=tmin,tmax=tmax,
                lpt=lpt,lpq=lpq,interstep=interstep):
    zmin = np.min(np.where(s.lpt/1e3>tmin))
    zmax = np.min(np.where(s.lpt/1e3>tmax))
    ts, te, F = parse_inter(lpt/1e3,lpq)
    plt.subplots()
    plt.plot(s.lpt/1e3,s.lpq/1e2)
    plt.plot(s.lpt[zmin:zmax]/1e3,s.lpq[zmin:zmax]/1e2)
    for i in range(np.size(F)):
        plt.plot([ts[i],te[i]],[F[i]/1e2,F[i]/1e2])


# define a function to save important results in hdf5 file
def sav(filename=filename, shot=shot, tcname=tcname, tau=tau, 
        dtau=dtau, Finter=Finter, Fmax=Fmax, efrac=efrac, A=A):
    overwrite = 0
    #open the file
    f = h5py.File(filename)
    #create subgroups
    tcloc = {9:'shelf',10:'floor',11:'floor'}
    if np.str(shot) not in f.keys():
        grp1 = f.create_group(np.str(shot))
        grp2 = grp1.create_group(tcloc[tcname])
    elif tcloc[tcname] not in f[np.str(shot)].keys():
        grp2 = f[np.str(shot)].create_group(tcloc[tcname])
    else:  # subgroups already exist and we want to overwrite
        overwrite = 1
        f[np.str(shot)][tcloc[tcname]]['elm_start[s]'][:] = tau
        f[np.str(shot)][tcloc[tcname]]['elm_duration[s]'][:] = dtau
        f[np.str(shot)][tcloc[tcname]]['elm_hf_only[W/m2]'][:] = \
                                        Fmax*efrac
        f[np.str(shot)][tcloc[tcname]]['hf_during_elm[W/m2]'][:] = \
                                        Fmax*efrac+Finter
        f[np.str(shot)][tcloc[tcname]]['energy_during_elm[J]'][:] = \
                                        (Fmax*efrac+Finter)*A*dtau
    #create data sets
    if overwrite == 0:
        grp2.create_dataset('elm_start[s]',data=tau)
        grp2.create_dataset('elm_duration[s]',data=dtau)
        grp2.create_dataset('elm_hf_only[W/m2]',data=Fmax*efrac)
        grp2.create_dataset('hf_during_elm[W/m2]',
                            data=Fmax*efrac+Finter)
        grp2.create_dataset('energy_during_elm[J]',
                            data=(Fmax*efrac+Finter)*A*dtau)
    #close the file
    f.close()

