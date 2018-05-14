# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:41:31 2017

@author: bartonjo

Take temperature history from TC data and calculate the heat flux
directly by assuming 1-D semi-infinite BCs, rearranging the temperature
history solution for the heat flux after taking the time derivative 
(F = dT/dt * ...), and using t = delta_t to get the heat flux after 
each time step.
If that works, then use this F to find the ELM heat flux using the 
techniques from hf_decon_v2.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import import_hfdecon_data as decondata
from cookbook import savitzky_golay as smooth
#from scipy.interpolate import interp1d
from int_erfc import i2erfc, i3erfc
import h5py
import pandas as pd

#---------------------------------------------------------------------
# User defined parameters
#

# hdf5 file name
filename = 'data.h5'

# data attributes from files
shot = 167277  #167408
tcname = 11 #10     # 9, 10, 11 for metal rings
pname = 17 #19      # 51, 19, (21), 15, (17)
tmin = 2.5     
tmax = 4. #3.9   # tmin and tmax is the time window for averaging data
de = 1. # tmin+de is time to start looking at ELM D-alpha data
cu = 0  # tmax-cu for ELM data to stop short of end-of-shot trans
mintime = 2.
maxtime = tmax  #3.9
septime = 1  # allowable time between elms; otherwise combine [ms]

# choose BCs
finite = 0      # 0 semi-infinite, 1 finite BCs

# Thickness of slab, d, and location of T analysis, x
d = .75e-2                      # measurement depth [m]
angle = 0                       # angle of button from horizontal
# .75e-2 insert TC, 2.77e-3 (3.18) 15o angled button, 1.57e-3 button
if angle>0:
    d += .75*6e-3*np.tan(np.radians(angle))
x = d                           # [m] from surface
d = 9.53e-3

# force the values of the heat flux if =1
force = 0
#Fmax = 1.24e6 #infinite 23.84e6  #finite 1.24
#Feff = 1.56e6 #infinite 3.16e6 #finite 1.56

# Material constants for pure Tungsten samples
#rho = 19.25e3                   # density [kg/m3]
#cp = 0.134e3                    # specific heat [J/kg/K]  0.134 J/g/K
#k = 125.                        # thermal conductivity [W/m/K]
                                # 173 RT, 125 (500C), 112 (927C)
# Material constants for Mo
rho = 10.2e3                     # density [kg/m3]
cp = .23e3                       # specific heat [J/kg/K]
k = 126.4                 # TZM thermal conductivity [W/m/K] at 100 C
# thermal diffusivity [m2/s]
kappa = k/rho/cp                

# other numerical parameters
num = range(20)                 # summation limit for finite BCs
t = np.linspace(1e-4,3.5,5000)     # seconds

#---------------------------------------------------------------------
# Extract useful parameters from the data
#

s = decondata.getdecondata(shot=shot, tcname=tcname, pname=pname)
 
# combine the false "double" ELMs into single ELMs
def elm_fix(s=s,d=10):
    estart = [s.elmstart[i] for i in range(np.size(s.elmstart)-1) 
            if s.elmstart[i] in s.dwt] #s.elmstart
    # find the duration points in the elmstart array that are in dwmhdf
    edur = [s.dur[i] for i in range(np.size(s.dur)-1) 
            if s.elmstart[i] in s.dwt]
#    edur = s.dur
    eend = [s.elmstart[i]+s.dur[i] for i in range(np.size(s.dur)-1) 
            if s.elmstart[i] in s.dwt]
    size = [s.dw[i] for i in range(np.size(s.dw)) 
            if s.dwt[i] in estart] #s.esize
    m = [i for i in range(np.size(eend)-1) if (estart[i+1]-eend[i])<d]
    m = np.array(m)
    elmstart = [estart[i] for i in range(np.size(eend)) \
            if i not in (m+1)]
    esize = [size[i] for i in range(np.size(eend)) if i not in (m+1)]
    dur1 = np.zeros(np.size(edur))
    for i in range(np.size(eend)-1):        
        if i not in m:
            dur1[i] = edur[i]
        else:
            dur1[i] = eend[i+1]-estart[i]
    dur = [dur1[i] for i in range(np.size(eend)) if i not in (m+1)]
    elmstart = np.array(elmstart)
    dur = np.array(dur)
    esize = np.array(esize)
#    plt.plot(s.elmstart,s.esize,'o-')
#    plt.plot(elmstart,esize,'o')
#    plt.plot(elmstart+dur,esize,'o')
    return elmstart, dur, esize

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

# Find the time offset when the strike point is on the TC
#tcloc = {9:1.41,10:1.354,11:1.324}
#osp = np.array([s.osp[i] for i in range(np.size(s.osp)) \
#            if s.osp[i]>0 and s.ospt[i]>2000.])
#ospt = [s.ospt[i] for i in range(np.size(s.osp)) \
#            if s.osp[i]>0 and s.ospt[i]>2000.]
#ospt = np.array(ospt)
#z1 = np.min(np.where(osp<tcloc[tcname+0]))
#z2 = np.min(np.where(osp[(z1+1):]<tcloc[tcname+0]))
#z3 = 1+z1+z2
#toffset = ospt[z3]/1e3  # convert ms to s
#toffset = toffset + t[0]   # simulation starts at t[0] not zero

# Find the wetted area of the ring
# average OSP radius
#z4 = np.max(np.where(ospt/1e3<tmax))
rloc = {9:1.38,10:1.3,11:1.3}
r = rloc[tcname]  #np.mean(osp[z3:z4+1])
print('mean radius = ', r)
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
    
# get the inter-ELM heat flux data
zmin = np.min(np.where(s.lpt/1e3>tmin+de))
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
try:
    lpt,lpq1,let,leq = elm_sep(s.lpt[zmin:zmax+1],s.lpq[zmin:zmax+1])
except:
    print('Cannot use the elm_sep function for this shot')
    lpt = s.lpt[zmin:zmax+1]
    lpq1 = s.lpq[zmin:zmax+1]
    let = 0.
    leq = 0.
lpt = np.array(lpt)
lpq1 = np.array(lpq1)
let = np.array(let)
leq = np.array(leq)
lpq = smooth(lpq1,151,3)
try:
    lptf,lpqf,letf,leqf = elm_sep(s.lpt,s.lpq)
except:
    print('Cannot use the elm_sep function for this shot')
    lptf = s.lpt
    lpqf = s.lpq
    letf = 0.
    leqf = 0.
lptf = np.array(lptf)
lpqf = np.array(lpqf)
letf = np.array(letf)
leqf = np.array(leqf)
#lpqf = smooth(lpqf,501,3)
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
Finter = lpqinter*1e6     # converted to W/m2
Fintermed = lpqavemed*1e6

# get ir average
zmin = np.min(np.where(s.irt/1e3>tmin+de))
zmax = np.min(np.where(s.irt/1e3>tmax))
try:
    irave = np.mean(s.irq2[zmin:zmax+1])/1e2  # MW/m2
except:
    irave = np.mean(s.irq[zmin:zmax+1])/1e2  # MW/m2
hirave = '{:0.2f}'.format(irave)
print('\nIR mean hf = '+hirave+' MW/m2 (raw average)')
#irt,irq,iet,ieq = elm_sep(s.irt[zmin:zmax+1],s.irq[zmin:zmax+1])
#irt = np.array(irt)
#irq = np.array(irq)
#irave = np.mean(irq)/1e2  # MW/m2
#hirave = '{:0.2f}'.format(irave)
#print('IR mean hf = '+hirave+' MW/m2 (ELM data removed)')

# for ease of plot manipulation and fitting, decrease the TC array size
tct = s.tct
tcTemp = s.tcTemp

#---------------------------------------------------------------------
# Define functions for model
#

# calculate the integral of the error function
def ierfc(z):
    return z*erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)

# calculate the heavyside function
def hside(c):
    return np.piecewise(c, [c>=0, c<0], [1,0])
    
# calculate the proportionality factor for F propto dT/dt
def get_factor(x=x,t0=tct[1]-tct[0],k=k,kappa=kappa,finite=finite,
               num=num,d=d):
    tstep = np.linspace(t0,2,1000)
    if finite>0:
        c1,c2 = 0,0
        for n in num:
            c1 += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*tstep)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*tstep))      
            c2 += erfc((2*d*(n+1)-x)/np.sqrt(4*kappa*tstep)) * \
                    ((2*d*(n+1)-x)/np.sqrt(16*kappa*tstep**3)) + \
                    erfc((2*d*n+x)/np.sqrt(4*kappa*tstep)) * \
                    ((2*d*n+x)/np.sqrt(16*kappa*tstep**3))
        dTdt_F = np.sqrt(kappa/k**2/tstep) * c1 + \
                np.sqrt(4*kappa*tstep/k**2) * c2
        dt = np.float(tstep[np.where(dTdt_F == np.max(dTdt_F))])
        c1,c2 = 0,0
        for n in num:
            c1 += -ierfc((2*d*(n+1)-x)/np.sqrt(4*kappa*dt)) - \
                    ierfc((2*d*n+x)/np.sqrt(4*kappa*dt))      
            c2 += erfc((2*d*(n+1)-x)/np.sqrt(4*kappa*dt)) * \
                    ((2*d*(n+1)-x)/np.sqrt(16*kappa*dt**3)) + \
                    erfc((2*d*n+x)/np.sqrt(4*kappa*dt)) * \
                    ((2*d*n+x)/np.sqrt(16*kappa*dt**3))
        dum = np.sqrt(kappa/k**2/dt) * c1 + \
                np.sqrt(4*kappa*dt/k**2) * c2
        factor = 1./dum
    else:
        # get the time step
        dTdt_F = np.sqrt(kappa/k**2/np.pi/tstep) * \
                np.exp(-x**2/4/kappa/tstep)
        dt = np.float(tstep[np.where(dTdt_F == np.max(dTdt_F))])
        factor = np.sqrt(k**2/kappa*np.pi*dt) * np.exp(x**2/4/kappa/dt)
    return factor, dt

#---------------------------------------------------------------------
# Find TC heat flux
#

# smooth the TC data
    # using 1500 points (.15 seconds) fit with 3rd order poly
That = smooth(s.tcTemp,1501,3) # TC10 is more noisey-use 5501 points
# get the numerical derivative of the TC Temp history
tct = tct/1e3
dTdt = np.gradient(That)/np.gradient(tct)
# get the TC heat flux measurement
factor, dt = get_factor()
F =  factor * dTdt 
Fhat = smooth(F,1501,3)

# get TC hf average
zmin = np.min(np.where(tct>tmin))
zmax = np.min(np.where(tct>tmax))
tcave = np.mean(F[zmin:zmax+1])/1e6  # MW/m2
htcave = '{:0.2f}'.format(tcave)
print('\nTC mean hf = '+htcave+' MW/m2')
#Finter = tcave*1e6

# get TC calculated ELM heat fluxes
zmin = np.min(np.where(tct>tmin+de))
zmax = np.min(np.where(tct>tmax))
q = (F[zmin:zmax+1]-Finter)/1e6
q1 = np.array([q[i] for i in range(np.size(q)) if q[i]>0])
t1 = np.array([tct[zmin+i] for i in range(np.size(q)) if q[i]>0])
dint = tct[1]-tct[0]
Q = dint * np.sum(q1)
    # get ELM properties
try:
    st, du, si = elm_fix(d=septime)
    zmin = np.min(np.where(st/1e3>tmin+de))
    ###try:
    zmax = np.max(np.where(st/1e3<tmax-cu))
    #zmin=np.min(np.where(s.elmstart/1e3>tmin+de))
    #zmax=np.max(np.where(s.elmstart/1e3<tmax-cu))
    ##except ValueError:
    ##    print('\nUsing ELM data to the end of the shot')
    ##    zmax = np.size(s.elmstart) - 1
        # find relative energy fraction of each ELM
    #esize = s.esize[zmin:zmax+1]
    esize = si[zmin:zmax+1]
    emax = np.max(esize)
    #dtau = s.dur[zmin:zmax+1]/1e3       # seconds
    dtau = du[zmin:zmax+1]/1e3       # seconds
    efrac = esize/emax
    #Fmax = Q/np.sum(dtau*efrac)*1e6  # for rectangle shape
    Fmax = Q/np.sum(.5*dtau*efrac)*1e6  # for triangle shape
    htcelm = '{:0.2f}'.format(Fmax/1e6)
    print('\nFmax = '+htcelm+' MW/m2')
        #get ELMs for whole shot
    #zmin = np.min(np.where(s.elmstart/1e3>tmin))
    zmin = np.min(np.where(s.elmstart/1e3>mintime))
    zmax = np.min(np.where(s.elmstart/1e3>maxtime)) 
        # ELM events
    #taustart = s.elmstart[zmin:zmax+1]/1e3   # seconds
    taustart = st[zmin:zmax+1]/1e3   # seconds
        # ELM duration
    #dtau = s.dur[zmin:zmax+1]/1e3       # seconds
    dtau = du[zmin:zmax+1]/1e3       # seconds
    tauend = taustart + dtau
    #esize = s.esize[zmin:zmax+1]
    esize = si[zmin:zmax+1]
    efrac = esize/emax
    Fdecon = Fintermed+np.zeros((np.size(tct)))
    for m in range(np.size(taustart)-1):
        felm = (hside(tct-taustart[m])-hside(tct-tauend[m])) * \
                efrac[m] * Fmax * \
                (1 - (tct-taustart[m])/(tauend[m]-taustart[m]))
        Fdecon += felm
    helmave = '{:0.2f}'.format(Fmax/1e6*np.mean(efrac))
    print('<Felm> = '+helmave+' MW/m2')
except:
    print('Cannot calculate ELM data for this shot')
    Fdecon = Fintermed+np.zeros((np.size(tct)))
    pass
## get the second derivative of the temp for the first order hf term
## take only the times when ELMs occur
#delay=.005
#for m in range(np.size(taustart)):
#    felm = (hside(tct+delay-taustart[m])-hside(tct+delay-tauend[m]))
#    if m == 0:
#        Fd = felm
#        delt = felm*(tct+delay-taustart[m])
#    else:
#        Fd += felm
#        delt += felm*(tct+delay-taustart[m])
#That2 = smooth(s.tcTemp,101,3)
#der = np.gradient(That2)/np.gradient(tct)
#d2Tdt2 = np.gradient(smooth(der,101,3))/np.gradient(tct)
## first order term has t-t0 factor....
##F2 = F + 541*.096*smooth(d2Tdt2,101,3)
##F2 = smooth(F2,101,3)  # examine trends lasting longer than 10 ms
#F2 = F + dTdt*Fd/4e-7 

#---------------------------------------------------------------------
# plotting
#

def compare_hfs(s=s,tct=tct,Fdecon=Fdecon,F=F,
                lpqf=lpqf,lptf=lptf):
    plt.figure()
    #plt.plot(tct,Fdecon/1e6,label='TC ELM deconvolution')
    #plt.plot(tct-.098,F2/1e6,'k', label='TC2 measurement')
    plt.plot(tct,F/1e6,'b', label='TC measurement')
    plt.plot(s.lptmed/1e3,s.lpqmed/1e2,'og', 
             label='LP (ELMs filtered out)')
    lps = smooth(s.lpq,151,3)
    plt.plot(s.lpt[0:3950]/1e3,lps[0:3950]/1e2,'g',
             label='LP (raw)', linewidth=2)
    #lpqfs = smooth(lpqf,501,3)
    #plt.plot(lptf/1e3,lpqfs/1e2)
    try:
        plt.plot(s.irt/1e3,s.irq2/1e2*1.4,'r',label='IRTV (peak)')
    except:
        print('plotting (r) integrated irtv data')
        plt.plot(s.irt/1e3,s.irq/1e2*1.4,'r',label='IRTV (int)')
    plt.xlim([0,6])
    plt.ylim([0,15])
    plt.legend(loc='best')
    plt.show()
    return

def compare_elms(s=s,tct=tct,Fdecon=Fdecon):
    plt.figure()
    plt.plot(tct,Fdecon/1e6)
    plt.plot(s.irt/1e3,s.irq2/1e2*1.4)
    plt.show()
    return

compare_hfs()
#compare_elms()


#---------------------------------------------------------------------
# saving data
#

# define a function to save important results in hdf5 file
def sav_old(filename=filename, shot=shot, tcname=tcname, tau=taustart, 
        dtau=dtau, Finter=Finter, Fmax=Fmax, efrac=efrac, A=A,
        tct=tct,Fdecon=Fdecon):
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
#    else:  # subgroups already exist and we want to overwrite
#        overwrite = 1
#        f[np.str(shot)][tcloc[tcname]]['elm_start[s]'][:] = tau
#        f[np.str(shot)][tcloc[tcname]]['elm_duration[s]'][:] = dtau
#        f[np.str(shot)][tcloc[tcname]]['elm_hf_only[W/m2]'][:] = \
#                                        Fmax*efrac
#        f[np.str(shot)][tcloc[tcname]]['hf_during_elm[W/m2]'][:] = \
#                                        Fmax*efrac+Finter
#        f[np.str(shot)][tcloc[tcname]]['energy_during_elm[J]'][:] = \
#                                        (Fmax*efrac+Finter)*A*dtau
    #create data sets
    if overwrite == 0:
        grp2.create_dataset('elm_start[s]',data=tau)
        grp2.create_dataset('elm_duration[s]',data=dtau)
        grp2.create_dataset('elm_hf_only[W m-2]',data=Fmax*efrac)
        grp2.create_dataset('hf_during_elm[W m-2]',
                            data=Fmax*efrac+Finter)
        grp2.create_dataset('energy_during_elm[J]',
                            data=(Fmax*efrac+Finter)*A*dtau)
        grp2.create_dataset('heat_flux_time[s]',data=tct)
        grp2.create_dataset('heat_flux_history[W m-2]',data=Fdecon)
    #close the file
    f.close()

def sav(filename=filename, shot=shot, tcname=tcname, tau=taustart, 
        dtau=dtau, Finter=Fintermed, Fmax=Fmax, efrac=efrac, A=A,
        tct=tct,Fdecon=Fdecon):
    # determine shelf or floor
    tcloc = {9:'shelf_',10:'floor_',11:'floor_'}
    # create key for file
    key = tcloc[tcname] + str(shot)
#    key = 'key2file'
    # create dataframe
    df1 = pd.DataFrame({'elm_start_s':tau, 'elm_duration_s':dtau,
            'elm_hf_only_Wm-2':Fmax*efrac,
            'hf_during_elm_Wm-2':Fmax*efrac+Finter,
            'energy_during_elm_J':(.5*Fmax*efrac+Finter)*A*dtau})
    df2 = pd.DataFrame({'heat_flux_time_s':tct,
            'heat_flux_history_Wm-2':Fdecon})
    df = pd.concat([df1,df2], axis=1)
    df.to_hdf(filename,key)
    return df
    





#-----------test code

# see notes for explanations
        
def fe(te,d=d,kappa=kappa):
    p1 = np.exp(-d**2/4/kappa/te)/np.sqrt(np.pi*te)
    p2 = d**2/4/kappa/te**2 - 1./2./te
    return p1*p2

def ge(te,taun,d=d,kappa=kappa):
    y = d/np.sqrt(2*kappa*te)
    p1 = 6./taun/np.sqrt(te)
    return p1*(i3erfc(y)-y*i2erfc(y)+y**2/3.*ierfc(y))


def firstd(ti,d=d,kappa=kappa,k=k):
    y = d/np.sqrt(4*kappa*ti)
    p1 = np.sqrt(kappa/k**2/np.pi/ti)
    return p1*np.exp(-y**2)
    
def secd(ti,d=d,kappa=kappa,k=k):
    p1 = d**2/4./kappa/ti**2 - 1./2./ti
    return firstd(ti,d=d,kappa=kappa,k=k)*p1

def coeffs(td,dur=2.5e-3):
    ts1 = td
    ts2 = td + dur
    te2 = td
    dTs1 = firstd(ts1)
    ddTs1 = secd(ts1)
    dTs2 = firstd(ts2)
    dTe2 = firstd(te2)
    ddTs2 = secd(ts2)
    ddTe2 = secd(te2)
    denom = dTs1*(ddTs2 - ddTe2) - ddTs1*(dTs2 - dTe2)
    A = (ddTs2 - ddTe2)/denom
    B = -(dTs2 - dTe2)/denom
    return A, B
    























