# -*- coding: utf-8 -*-
"""
Finite slab and semi-infinite 1-D thermal model for an average heat
flux and a combination of inter-ELM and ELM heat fluxes on the surface.

Created on Thu Jun 22 12:59:41 2017

@author: bartonjo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc, erf
#import import_hfdecon_data as decondata
#from int_erfc import i3erfc

#---------------------------------------------------------------------
# Material constants
## pure Tungsten samples
#rho = 19.25e3                   # density [kg/m3]
#cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
#k = 125.                        # thermal conductivity [W/m/K]
                                # 173 RT, 125 (500C), 112 (927C)
## SiC
#rho = 3.21e3                    # density [kg/m3]
#cp = .64e3                      # specific heat [J/kg/K]  
#k = 300.                        # thermal conductivity [W/m/K]

# ATJ graphite (@ 100 C)
rho = 1758.9                    # density [kg/m3]
cp = 1125.3                      # specific heat [J/kg/K]  
k = 87.19                        # thermal conductivity [W/m/K]

kappa = k/rho/cp                # thermal diffusivity [m2/s]

# Thickness of slab, d, and location of T analysis, x
d = 2e-3#6e-3                      # depth [m]
angle = 0                       # 1 = 15 deg angle button depth
if angle>0:
    d += 1.2e-3 #.75*6e-3*np.tan(np.radians(angle))
x = d
# TC tolerance temperatures
if x>0:
    if angle>0:
        tol = 982
    else:
        tol = 593
else:
    tol = 0

#---------------------------------------------------------------------
# Simulation constants
# heat fluxes arbitrary shot
Finter = 1.5e6 #10e6            # inter-ELM [W/m2]
Fmax = 8e6 #5.*Finter           # max ELM 
Feff = 5e5 #2.19e6           # average heat flux
## heat fluxes for worst case 167297
#Finter = 17e6            # inter-ELM [W/m2]
#Fmax = 39e6           # max ELM 
#Feff = 19e6           # average heat flux
#print('parallel heat flux = ', Feff/np.sin(np.radians(2.1))/1e6)

# choose BCs
finite = 1      # 0 semi-infinite, 1 finite BCs

# time parameters
tend = 4.           # [s]
tstart = 1e-4
tpoints = 55000

## get ELM start and end times as well as their relative sizes
## or synthetically imput these values
## from data:
#shot = 167353
#tcname = 9
#pname = 51
#s = decondata.getdecondata(shot=shot, tcname=tcname, pname=pname)
s=1.
tcname=1
# synthetic: (Assume all ELMs are the same size = Fmax)
freq = 35 #10.           # ELM frequency [Hz]
dur = 2.5e-3 #5e-3          # ELM duration [s]
ts = np.arange(tstart+1./freq,tend,1./freq)
te = ts + dur
efrac = np.ones(np.size(te))


#---------------------------------------------------------------------
# Initialize arrays
num = range(20)                 # summation limit for finite BCs
t = np.linspace(tstart,tend,tpoints)     # seconds
prof = np.linspace(0,d,100)      # depth profile [m]
tprof = np.linspace(tstart,tend,np.int(tpoints/100.))

#---------------------------------------------------------------------
# Thermal model functions
# calculate the integral of the error function
def ierfc(z):
    return z*erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)

# calculate the heavyside function
def hside(c):
    return np.piecewise(c, [c>=0, c<0], [1,0])
    
# calculate the temperature due to a single constant heat flux
def T_chf(F,x=x,t=t,k=k,kappa=kappa,finite=finite,num=num,d=d):
    if finite>0:
        c = 0
        for n in num:
            c += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*t)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*t))      
        T = F/k * np.sqrt(4*kappa*t) * c
    else:
        T = F/k * np.sqrt(4*kappa*t) * (-ierfc(x/np.sqrt(4*kappa*t)))
    return T

# calculate the temperature due to a single constant heat flux + ELMs
def T_ehf(Finter,Fmax,x=x,t=t,k=k,kappa=kappa,finite=finite,num=num,
          d=d,ts=ts,te=te,efrac=efrac):
    T = T_chf(Finter,x=x,t=t,k=k,kappa=kappa,finite=finite,
                  num=num,d=d)
    if finite>0:
        for m in range(np.size(ts)):
            dts = t-ts[m]
            dte = t-te[m]
            dts = np.array([i if i>0 else 1 for i in dts])
            dte = np.array([i if i>0 else 1 for i in dte])
            # triangle solution
            #...
            # rectangle solution
            c1,c2 = 0,0
            for n in num:
                c1 += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*dts)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*dts))
                c2 += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*dte)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*dte))
            telm = Fmax/k * efrac[m] * \
                (hside(t-ts[m]) * np.sqrt(4*kappa*dts) * c1 -
                hside(t-te[m]) * np.sqrt(4*kappa*dte) * c2)
            T += telm
    else:
        for m in range(np.size(ts)):
            dts = t-ts[m]
            dte = t-te[m]
            te2 = ts[m]+.5e-3
            dte2 = t-te2
            dts = np.array([i if i>0 else 1 for i in dts])
            dte = np.array([i if i>0 else 1 for i in dte])
            dte2 = np.array([i if i>0 else 1 for i in dte2])
            telm = np.zeros(np.size(t)) 
            # solution for electrons?
#                # delta function
##            telm += Fmax/5./50./k * efrac[m] * hside(t-ts[m]) * \
##                np.sqrt(kappa/np.pi/dts) * np.exp(-x**2/4./kappa/dts)
#                # rectangle
#            telm += Fmax*20./k * efrac[m] * \
#                (-hside(t-ts[m]) * np.sqrt(4*kappa*(dts)) * 
#                ierfc(x/np.sqrt(4*kappa*(dts))) +
#                hside(t-te2) * np.sqrt(4*kappa*(dte2)) * 
#                ierfc(x/np.sqrt(4*kappa*(dte2))))
            #triange solution
#            b = Fmax/k * efrac[m] * 1./(te[m]-ts[m])
#            telm += Fmax/k * efrac[m] * hside(t-ts[m]) * \
#                np.sqrt(4*kappa*dts) * (-ierfc(x/np.sqrt(4*kappa*dts)))
#            telm += -b * hside(t-ts[m]) * \
#                np.sqrt(kappa) * (4*dts)**1.5 * \
#                (-i3erfc(x/np.sqrt(4*kappa*dts)))
#            telm += b * hside(t-te[m]) * \
#                np.sqrt(kappa) * (4*dte)**1.5 * \
#                (-i3erfc(x/np.sqrt(4*kappa*dte)))
            # rectangle solution
            telm += Fmax/k * efrac[m] * \
                (-hside(t-ts[m]) * np.sqrt(4*kappa*(dts)) * 
                ierfc(x/np.sqrt(4*kappa*(dts))) +
                hside(t-te[m]) * np.sqrt(4*kappa*(dte)) * 
                ierfc(x/np.sqrt(4*kappa*(dte))))
            T += telm
    return T

# return the ELM/inter-ELM heat flux BC time history
def hf_elms(Finter,Fmax,t=t,ts=ts,te=te,efrac=efrac):
    F = np.zeros(np.size(t)) + Finter
    for m in range(np.size(ts)):
        felm = np.zeros(np.size(t))
        # triangle elms
#        felm += (hside(t-ts[m])-hside(t-te[m])) * \
#                efrac[m] * Fmax * (1 - (t-ts[m])/(te[m]-ts[m]))
        # rectangle elms
        felm += (hside(t-ts[m])-hside(t-te[m])) * \
                efrac[m] * Fmax
        # electron contribution?
#        te2 = ts[m] + .5e-3
#        felm += (hside(t-ts[m])-hside(t-te2)) * \
#                efrac[m] * Fmax *20
        F += felm
    return F

# return the temperature profile given an initial profile and no heat
# flowing from the ends
def relax(m,b,prof=prof,t=tprof,k=k,kappa=kappa,num=num,d=d):
    # m is the slope and is assumed to be positive, b is T_surface: 
    # Tinit = -m*x + b
    T = np.zeros([np.size(prof),np.size(t)])    
    for x in prof:
        s1,s0 = np.zeros(np.size(t)),np.zeros(np.size(t))
        for n in num:
            if n>0:
                s1 += ierfc((2*n*d-x)/np.sqrt(4*kappa*t))
            s0 += -ierfc(((2*n+1)*d-x)/np.sqrt(4*kappa*t))
            s0 += ierfc((2*n*d+x)/np.sqrt(4*kappa*t))
            s0 += -ierfc(((2*n+1)*d+x)/np.sqrt(4*kappa*t))
        temp = m*np.sqrt(4*kappa*t)*(s1+s0) + (-m*x+b)
        T[np.where(prof==x)[0][0]][:] = temp
    return T

# return the temperature profiles of a composite material (x into ATJ),
# where the temperatures of the two materials start out uniformly const
def comp_relax(Tw,Tg,H=1.,x=prof,tprof=tprof,k=k,kappa=kappa):
    # initial temperature of x material to reference ATJ at 0 temp   
    Ti = Tw-Tg  
    # graphite (ATJ) properties at 100 C
    kg = 87.19  # [W/m/K]
    rhog = 1758.92  #[kg/m3]
    cpg = 1125.30  # [J/kg/K]
    kappag = kg/rhog/cpg
    coeff = k/np.sqrt(kappa)*Ti/(k/np.sqrt(kappa) + kg/np.sqrt(kappag))
    T = np.zeros([2*np.size(x),np.size(tprof)])
    if H == np.inf:
        for t in tprof:
            T1 = coeff * (1 + kg/np.sqrt(kappag)/(k/np.sqrt(kappa)) *
                    erf(x/np.sqrt(4*kappa*t)))
            T2 = coeff * erfc(x/np.sqrt(4*kappag*t))
            T[:,np.where(tprof==t)[0][0]] = np.append(np.flipud(T2),T1)
    else:
        h = H*(k/np.sqrt(kappa)+kg/np.sqrt(kappag))/(k*kg/
                np.sqrt(kappag))
        hg = h * (k*kg/np.sqrt(kappag))/(kg*k/np.sqrt(kappa))
        for t in tprof:
            T1 = coeff*(1 + kg/np.sqrt(kappag)/(k/np.sqrt(kappa))*
                    (erf(x/np.sqrt(4*kappa*t)) + 
                    np.exp(h*x+h**2*kappa*t)*
                    erfc(x/np.sqrt(4*kappa*t)+h*np.sqrt(kappa*t))))
            T2 = coeff*(erfc(x/np.sqrt(4*kappag*t)) - 
                    np.exp(hg*x+hg**2*kappag*t)*
                    erfc(x/np.sqrt(4*kappag*t)+hg*np.sqrt(kappag*t)))
            T[:,np.where(tprof==t)[0][0]] = np.append(np.flipud(T2),T1)
    newx = np.append(-np.flipud(x),x)
    return newx, T

# return the distance from the measurement to the OSP
def OSPdist(x,tmin,tmax,s=s, tcname=tcname):
    # tmin, tmax in seconds. lengths in m.
    tclocs = {9:141., 10:135.4, 11:132.4}
    tc = tclocs[tcname]/100.
    zmin = np.min(np.where(s.ospt/1e3 > tmin))
    zmax = np.min(np.where(s.ospt/1e3 > tmax))
    y = np.array([np.abs(tc-osp) for osp in s.osp[zmin:zmax+1]])
    dist = np.sqrt(x**2+y**2)
    return dist, s.ospt[zmin:zmax+1]/1e3

#---------------------------------------------------------------------
# Plotting
## plot heat flux BC time history
## Finter=1e6,Fmax=1.2e6, dur =3ms, freq = 100
#plt.subplots()
#plt.plot(t,hf_elms(Finter,Fmax)/1e6, linewidth=2)
##plt.xlim([.99,1.05])
##plt.ylim([.95,2.25])
#plt.xlabel('Time [s]', fontsize=24)
#plt.ylabel(r'Heat flux [MW/m$^2$]', fontsize=24)
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.subplots_adjust(bottom=.12)
#plt.show()

## plot T history for both heat flux BCs
## finite BC: ang=1, Finter=4.8e6, Feff=5e6 dur=2.5ms, freq = 10
## infinte BC: ang=1, Finter=, Feff=, dur=2.5ms, freq=10
##Finter = 11.5e6*fact            # inter-ELM [W/m2]
##Fmax = 1.7*Finter           # max ELM 
##Feff = 12e6*fact 
#plt.subplots()
##plt.plot(t,T_ehf(Finter,Fmax)+100, label=r'$q_{inter}$ and $q_{ELM}$')
#plt.plot(t,T_chf(Feff), label = 'Finite BCs')
#plt.plot(t,T_chf(Feff, finite=0), label = 'Semi-infinite BCs')
##plt.plot([0,4],[3422,3422],'k')
#plt.plot([2.5,2.5],[0,200],'k')
##plt.plot([0,4],[1000,1000],'--k')
##plt.plot([0,4],[tol,tol],'r')
##plt.xlim([.99,1.05])
##plt.ylim([.95,2.25])
#plt.xlabel('Time [s]', fontsize=24)
#plt.ylabel(r'$\Delta$T', fontsize=24)
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.subplots_adjust(bottom=.12, left=.14)
#plt.legend(loc='best', fontsize=18)
#plt.show()

## plot the depth profile
## infinte BC: ang=1, Finter=, Feff=, dur=2.5ms, freq=10
#plt.subplots()
##plt.plot(t,T_ehf(Finter,Fmax,x=prof,t=3.5)+100)
#plt.plot(prof*1e3,T_chf(Feff,x=prof,t=2.0)+0,linewidth=2)
##plt.plot([0,d*1e3],[1000,1000],'--k')
#plt.plot([0,d*1e3],[tol,tol],'r')
##plt.xlim([.99,1.05])
##plt.ylim([800,1500])
#plt.xlim([0,7.2])
#plt.ylim([500,1100])
#plt.xlabel('Depth [mm]', fontsize=24)
#plt.ylabel(r'Temperature [C]', fontsize=24)
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.subplots_adjust(bottom=.12, left=.14)
#plt.show()

## plot the relaxation of the depth profile
## dT=570, T_surface=1100, d:angled sample
##Tr = relax((570.)/d,1100)
#plt.subplots()
#for i in range(6):
#    s = i*10
#    h = '{:0.2f}'.format(t[s])
#    plt.plot(prof*1e3,Tr[:,s], linewidth=2, label= h + ' sec')
#plt.xlim([0,7.2])
#plt.xlabel('Depth [mm]', fontsize=24)
#plt.ylabel(r'Temperature [C]', fontsize=24)
#plt.title('T relaxation after 2 sec shot \non angled sample', 
#          fontsize=24)
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.subplots_adjust(bottom=.12, left=.14, top=.88)
#plt.legend(loc='best', fontsize=18)
#plt.show()    

##plot the relaxation of T when W and graphite are in contact and have
## different temperatures. H=np.inf is perfect conductance b/w materials
#Tw = 900.
#Tg = 400.
#cond = np.inf  # conductivity b/w the two materials
#newx = np.linspace(0,.01,1000)  # 1 cm into each material
#newt = np.linspace(1e-4,5*60,550)  # 5 minutes
#rho = 19.25e3
#cp = .134e3
#k = 125.
#kappa = k/rho/cp
#x, T = comp_relax(Tw,Tg,H=cond,tprof=newt,x=newx,k=k,kappa=kappa)
#plt.figure()
#for i in range(6):
#    n = i*10
#    h = '{:0.2f}'.format(newt[n])
#    plt.plot(x*1e2,T[:,n]+Tg,linewidth=2, label=h + ' sec')
#plt.ylim([Tg-50,Tw+50])
#plt.xlabel('Depth [cm]', fontsize=24)
#plt.ylabel(r'Temperature [C]', fontsize=24)
#plt.title(
#    'T relaxation in W (x>0) and ATJ (x<0) \n(perfect conductance)',
#    fontsize=24)
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.subplots_adjust(bottom=.12, left=.14, top=.88)
#plt.legend(loc='best', fontsize=18)
#plt.show() 


##plot the 1-d predictions on top of the data from shot 172370
#import pandas as pd
#ir = pd.read_table('ir60_172370.dat', sep='\s+', names=['irt','irq'])
#irt = ir.irt
#irq = ir.irq
#tc = pd.read_table('tctemp04_172370.dat', sep='\s+', 
#                   names=['tct','tcTemp'])
#datt = tc.tct/1e3
#datTemp = tc.tcTemp
#f, (ax1,ax2) = plt.subplots(2, sharex=True)
#
#ax1.plot(irt,irq,'r')
#ax1.plot([2.05,2.05],[0,1000],'--k')
#ax1.plot([4.05,4.05],[0,1000],'--k')
#ax1.text(4.8,5.5,'inter-ELM ~ 1.5 MW/m2')
#ax1.text(4.8,5.1,'ELM ~ 4 MW/m2 @ 30 Hz')
#ax1.set_ylim([0, 6]) 
#ax1.set_ylabel(r'q$_\perp$ [MW/m$^2$]', fontsize=24)
#ax1.tick_params(axis='y', labelsize=18)
#
#ax2.plot(t+2.05,T_ehf(1.5e6*8, 4e6*8, x=d)+350, linewidth=2, 
#         label='model at 7.2 mm depth')
#ax2.plot(datt,datTemp, label='TC under middle button, 172370') 
#ax2.plot([2.05,2.05],[0,1000],'--k')
#ax2.text(2.05,900,'H-mode starts', fontsize=16)
#ax2.plot([4.05,4.05],[0,1000],'--k')
#ax2.text(4,500,'OSP moved off DiMES', fontsize=16)
#ax2.text(4.5,800,'TC saturated', fontsize=14, color='green') 
#ax2.set_xlim([1, 6])
#ax2.set_ylim([0, 1000])          
#ax2.set_xlabel('t [s]', fontsize=24)
#ax2.set_ylabel(r'$\Delta$T', fontsize=24)
#ax2.tick_params(axis='x', labelsize=18)
#ax2.tick_params(axis='y', labelsize=18)
#f.subplots_adjust(bottom=.12, left=.14, hspace=0.1)
#ax2.legend(loc='best', fontsize=18)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#plt.show()