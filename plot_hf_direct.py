# -*- coding: utf-8 -*-
"""
Created on Thu May  4 17:05:48 2017

@author: bartonjo

Plot data generated by the hf_direct.py script
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import import_hfdecon_data as decondata
from cookbook import savitzky_golay as smooth
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

#---------------------------------------------------------------------
# Compare IR data to TC deconvolution
# shot 167353, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5, de=1Elm
    # get inter-ELM IR data
irti, irqi, irte, irqe = elm_sep(s.irt,s.irq2)
    # smooth inter-ELM IR data
irqsm = smooth(irqi,1001,3)
    # interpolate smoothed data on full IR data scale
fint = interp1d(irti,irqsm)
irback = fint(s.irt)
    # subtract IR heat flux background to get ELM heat fluxes
irelms = (s.irq2 - irback)*1.
    # plot TC ELM heat flux and IR ELM heat flux
f, ax = plt.subplots()
ax.plot(s.irt/1e3,irelms/1e2, 'r', linewidth=2, label='IR ELMs')
ax.plot(tct,(Fdecon-Fintermed)/1e6, 'b', linewidth=2,
         label='TC ELM calculation')
ax.set_xlabel('t [s]', fontsize=24)
ax.set_ylabel(r'ELM heat flux [MW/m$^2$]', fontsize=24)
#plt.xlim([2.1,3.1])
#plt.xlim([2.6,3.1])
#plt.xlim([2.67,2.75])
plt.xlim([2.68,3.11])
ax.set_ylim([-.1,6.1])
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)
plt.legend(loc='upper left', fontsize=18)
plt.subplots_adjust(bottom=.12)
#### This is for the zoomed in inlay
###plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
###plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
###plt.xlim([2.68,3.11])
###ax.set_ylim([-.1,6.1])
##plt.show()

##---------------------------------------------------------------------
## T vs t from TC data before and after smoothing
## shot 167353
#plt.subplots()
#plt.plot(s.tct/1e3,s.tcTemp, linewidth=2, label='raw data')
#plt.plot(tct,That,'r',linewidth=1.5, label='smoothed data')
#for i in range(6):
#    plt.plot([2+.5*i,2+.5*i],[0,1000],'k',linewidth=.5)
#plt.xlabel('t [s]', fontsize=24)
#plt.ylabel('Temperature [C]', fontsize=24)
#plt.xlim([1.6,5])
#plt.ylim([40,220])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.legend(loc='upper left', fontsize=18)
#plt.subplots_adjust(bottom=.12)
#plt.show()

#---------------------------------------------------------------------
# heat flux vs t comparing LP and TC data for high and low ELM freqs
## Low ELM frequency: shot 167353, finite=0, Finter=lpqinter, shelf
#plt.subplots()
#plt.plot(tct-.0,F/1e6, label='TC', linewidth=1)
##lpqfs = smooth(lpqf,501,3)
##plt.plot(lptf[0:3950]/1e3,lpqfs[0:3950]/1e2,'r')
#lps = smooth(s.lpq,151,3)
#plt.plot(s.lpt[0:3950]/1e3,lps[0:3950]/1e2,'k',
#         label='LP', linewidth=2)
##plt.plot(lpt/1e3,lpq/1e2,
##         label='LP', linewidth=2)
##plt.plot(s.lptmed/1e3,s.lpqmed/1e2,'o', label='LP (inter-ELM)')
#plt.plot(s.irt/1e3,s.irq2/1e2, 'r', 
#         label='IRTV (peak value)', linewidth=2)
#plt.xlabel('t [s]', fontsize=24)
#plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=24)
#plt.xlim([1.8,5.4])
#plt.ylim([-.1,6.1])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
##plt.legend(loc='upper right', fontsize=18)
#plt.legend(bbox_to_anchor=(1.1, 1.1), fontsize=16)
#plt.subplots_adjust(bottom=.12)
#plt.show()
## High ELM frequency: shot 167320, finite=0, Finter=lpqinter, shelf
#plt.subplots()
#plt.plot(tct-.0,F/1e6, label='TC', linewidth=2)
#lps = smooth(s.lpq,151,3)
#plt.plot(s.lpt[0:3950]/1e3,lps[0:3950]/1e2,'k',
#         label='LP', linewidth=3)
##plt.plot(s.lptmed/1e3,s.lpqmed/1e2,'o', label='LP (ELM filtered)')
##plt.plot(s.irt/1e3,s.irq/1e2, 'm', 
##         label='IRTV (integrated across insert)')
#plt.plot(s.irt/1e3,s.irq2/1e2, 'r', 
#         label='IRTV (peak value)', linewidth=1)
#plt.xlabel('t [s]', fontsize=24)
#plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=24)
#plt.xlim([1.8,4.4])
#plt.ylim([-.1,8.1])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.legend(bbox_to_anchor=(1.1, 1.1), fontsize=16)
#plt.subplots_adjust(bottom=.12)
#plt.show()


##---------------------------------------------------------------------
## Illustration of ELM energy calculation to get ELM heat flux
## Stacked plots of F-Finter and F_ELMs vs t
## shot 167353, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5
#f, (ax1, ax2) = plt.subplots(2, sharex=True) 
#ax1.plot(tct, (F-Finter)/1e6)
#ax2.plot(tct,(Fdecon-Finter)/1e6, linewidth=2,
#         label='TC ELM calculation')
#ax1.set_ylim([0,2])
#ax2.set_ylim([0,4])
#ax2.set_xlim([2.3,3.3])
#ax1.set_ylabel(r'$q_{\perp} - q_{inter}$' '\n' '[MW/m$^2$]',fontsize=24)
#ax1.tick_params(axis='y', labelsize=18)
#ax2.set_ylabel(r'$q_{ELMs}$' '\n' ' [MW/m$^2$]',fontsize=24)
#ax2.tick_params(axis='y', labelsize=18)
#ax2.tick_params(axis='x', labelsize=18)
#ax2.set_xlabel('t [s]',fontsize=24)
## Fine-tune figure; make subplots close to each other and 
## hide x ticks for all but bottom plot.
#f.subplots_adjust(hspace=0, bottom=.12, left=.17)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)


##---------------------------------------------------------------------
## Sequence of TC shots that show bowing in the temp histories
## shots: 167431,32,39,40, need OSP, TC10, and TC11 data (no modeling)
#shots = [167431, 167432, 167439, 167440]
#i=0
#for shot in shots:
#    tcname = 10
#    pname = 19 
#    s = decondata.getdecondata(shot=shot, tcname=tcname, pname=pname)
#    if i==0:
#        ospt = np.zeros([np.size(shots),np.size(s.ospt)])
#        osp = np.zeros([np.size(shots),np.size(s.osp)])
#        tct10 = np.zeros([np.size(shots),np.size(s.tct)])
#        tcTemp10 = np.zeros([np.size(shots),np.size(s.tcTemp)])
#    try:
#        ospt[i,:] = s.ospt/1e3
#        osp[i,:] = s.osp*100
#    except ValueError:
#        ospt[i,:] = s.ospt[np.abs(np.size(s.ospt)-\
#                np.size(ospt[0,:])):]/1e3
#        osp[i,:] = s.osp[np.abs(np.size(s.osp)-np.size(osp[0,:])):]*100
#    try:
#        tct10[i,:] = s.tct/1e3
#        tcTemp10[i,:] = s.tcTemp
#    except ValueError:
#        print(i)
#        tct10[i,:] = s.tct[np.abs(np.size(s.tct)-np.size(tct10[0,:])):]/1e3
#        tcTemp10[i,:] = s.tcTemp10[np.abs(np.size(s.tcTemp)-\
#                np.size(tcTemp10[0,:])):]
#    tcname = 11
#    s = decondata.getdecondata(shot=shot, tcname=tcname, pname=pname)
#    if i==0:
#        tct11 = np.zeros([np.size(shots),np.size(s.tct)])
#        tcTemp11 = np.zeros([np.size(shots),np.size(s.tcTemp)])
#    try:
#        tct11[i,:] = s.tct/1e3
#        tcTemp11[i,:] = s.tcTemp
#    except ValueError:
#        tct11[i,:] = s.tct[np.abs(np.size(s.tct)-\
#                np.size(tct11[0,:])):]/1e3
#        tcTemp11[i,:] = s.tcTemp10[np.abs(np.size(s.tcTemp)-\
#                np.size(tcTemp11[0,:])):]
#    i+=1
## get nonbowing shot data
#s1 = decondata.getdecondata(shot=167408, tcname=10, pname=pname)
#s2 = decondata.getdecondata(shot=167408, tcname=11, pname=pname)
#f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True) 
#ax1.plot(s1.ospt/1e3, s1.osp*100, '--', linewidth=2, 
#         label='167408 (typical response)')
##ax1.plot(ospt[0,:], osp[0,:], linewidth=2, label='167431')
##ax1.plot(ospt[1,:], osp[1,:], linewidth=2, label='167432')
#ax1.plot(ospt[2,:], osp[2,:], linewidth=2, label='167439 (bowing)')
#ax1.plot(ospt[3,:], osp[3,:], linewidth=2, label='167440 (bowing)')
#st = 100
#ax2.plot(s1.tct[::st*50]/1e3,s1.tcTemp[::st*50],'--', 
#         linewidth=2)  # TC10 data
##ax2.plot(tct10[0,::st],tcTemp10[0,::st], linewidth=2)  # TC10 data
##ax2.plot(tct10[1,::st],tcTemp10[1,::st], linewidth=2)
#ax2.plot(tct10[2,::st],tcTemp10[2,::st], linewidth=2)
#ax2.plot(tct10[3,::st],tcTemp10[3,::st], linewidth=2)
#ax3.plot(s2.tct[::st]/1e3,s2.tcTemp[::st],'--', 
#         linewidth=2)  # TC11 data
##ax3.plot(tct11[0,::st],tcTemp11[0,::st], linewidth=2)  # TC11 data
##ax3.plot(tct11[1,::st],tcTemp11[1,::st], linewidth=2)
#ax3.plot(tct11[2,::st],tcTemp11[2,::st], linewidth=2)
#ax3.plot(tct11[3,::st],tcTemp11[3,::st], linewidth=2)
#ep = .25  #line thickness compensation
#ax1.plot([.8,.8],[132+ep,137-ep],color='0.50', linewidth=13)
#ax1.plot([0,7],[135.4,135.4],'k',linewidth=.5)
#ax1.plot([0,7],[132.4,132.4],'k',linewidth=.5)
#ax1.plot([5.65,5.65],[0,1000],'g',linewidth=.75)
#ax1.plot([5.65,5.65],[0,1000],'--r',linewidth=.75)
#ax2.plot([5.65,5.65],[0,1000],'g',linewidth=.75)
#ax2.plot([5.65,5.65],[0,1000],'--r',linewidth=.75)
#ax3.plot([5.65,5.65],[0,1000],'g',linewidth=.75)
#ax3.plot([5.65,5.65],[0,1000],'--r',linewidth=.75)
#ax1.plot([5.15,5.15],[0,1000],'b',linewidth=.75)
#ax2.plot([5.15,5.15],[0,1000],'b',linewidth=.75)
#ax3.plot([5.15,5.15],[0,1000],'b',linewidth=.75)
##ax2.plot([5.65,5.65],[0,1000],'k',linewidth=.5)
##ax3.plot([5.65,5.65],[0,1000],'k',linewidth=.5)
#ax1.set_ylim([130.4,137.4])
#ax2.set_ylim([60,395])
#ax3.set_ylim([60,310]) #510])
#ax3.set_xlim([.8,15.2])
#ax1.set_ylabel('OSP [cm]',fontsize=24)
#ax1.tick_params(axis='y', labelsize=18)
#ax2.set_ylabel('T [C]',fontsize=24)
#ax2.tick_params(axis='y', labelsize=18)
#ax3.set_ylabel('T [C]',fontsize=24)
#ax3.tick_params(axis='y', labelsize=18)
#ax3.tick_params(axis='x', labelsize=18)
#ax3.set_xlabel('t [s]',fontsize=24)
#ax1.legend(loc='upper right', fontsize=18)
## Fine-tune figure; make subplots close to each other and 
## hide x ticks for all but bottom plot.
#f.subplots_adjust(hspace=0, bottom=.12)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)


##---------------------------------------------------------------------
## Compare IR data to TC deconvolution
## shot 167353, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5, de=1?
### shot 167354, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5
## shot 167355, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5, de=1?
## shot 167356, finite=0, Finter=lpqinter, tmin=2.1, tmax=4.5, de=1?
#    # get inter-ELM IR data
#irti, irqi, irte, irqe = elm_sep(s.irt,s.irq2) # USE PEAK?
#    # smooth inter-ELM IR data
#irqsm = smooth(irqi,1001,3)
#    # interpolate smoothed data on full IR data scale
#fint = interp1d(irti,irqsm)
#irback = fint(s.irt)
#    # subtract IR heat flux background to get ELM heat fluxes
#irelms = s.irq - irback
#    # get ELM energies
#ire = np.zeros(np.size(taustart)-1)
#tce = np.zeros(np.size(taustart)-1)
#for m in range(np.size(taustart)-1):
#    # Integrate the IR data during ELMs to get energy of each ELM
#    zmin = np.min(np.where(s.irt/1e3>taustart[m]))
#    zmax = np.min(np.where(s.irt/1e3>tauend[m]))
#    dt = (s.irt[zmin+1]-s.irt[zmin])/1e3
#        # MJ (.45 m2 shelf ring area)
#    ire[m] = dt * np.sum(irelms[zmin:zmax+1])/1e2*.45   
#    # Integrate the TC model data to get the energy of each ELM
#        # MJ (.45 m2 shelf ring area)
#    tce[m] = dtau[m] * 1. *(Fmax*efrac[m])/1e6 * .45   
#
#    # remove zero data
#tce = np.array([tce[i] for i in range(np.size(ire)) if ire[i]>0])
#ire = np.array([ire[i] for i in range(np.size(ire)) if ire[i]>0])
##    # After set-up combine shots
##ire = np.concatenate((ire353,ire355,ire356))
##tce = np.concatenate((tce353,tce355,tce356))
#    # fit a line to the data
#def line(indep,m):
#    return m*indep
#popt, pcov = curve_fit(line,tce*1e3,ire*1e3)
#rsq = pearsonr(ire*1e3,line(tce*1e3,*popt))[0]
#print('Rsq = ', rsq)
#print('m = ', np.float(popt))
#    # plot TC ELM energy and IR ELM energy
#plt.subplots()
#plt.errorbar(tce*1e3,ire*1e3, xerr=.7, yerr=.2, fmt='+')
##plt.plot(tce*1e3,ire*1e3,'o')
#plt.plot(np.linspace(0,20),line(np.linspace(0,20),*popt),
#         linewidth=3)
##plt.plot(np.linspace(0,20),line(np.linspace(0,20),1),'k')
#plt.xlabel('TC ELM energy [kJ]', fontsize=24)
#plt.ylabel('IR ELM energy [kJ]', fontsize=24)
#plt.xlim([0,20])
#plt.ylim([0,20])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
##plt.legend(loc='upper right', fontsize=18)
#plt.subplots_adjust(bottom=.12)
#plt.show()

##---------------------------------------------------------------------
## Show deconvolution results
## shot 167432, tcname = 11, pname = 17, tmin = 2.1, tmax = 4.5 
#
#plt.subplots()
#plt.plot(tct,Fdecon/1e6,'c',label='TC ELM deconvolution', 
#         linewidth=1.5)
#plt.plot(tct,F/1e6,'b', label='TC measurement', linewidth=2)
#plt.plot(s.lptmed/1e3,s.lpqmed/1e2,'og', 
#         label='LP (ELMs filtered out)')
#plt.xlabel('t [s]', fontsize=24)
#plt.ylabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=24)
#plt.xlim([3.4,4.5])
#plt.ylim([1.5,11.5])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.legend(loc='upper left', fontsize=18)
#plt.subplots_adjust(bottom=.12)
#plt.show()

##---------------------------------------------------------------------
## plot how hf affects surface and leading edge temperature
## Material constants for TZM
#rho = 10.2e3                     # density [kg/m3]
#cp = .23e3                       # specific heat [J/kg/K]
#k = 128.                         # thermal conductivity [W/m/K]
## thermal diffusivity [m2/s]
#kappa = k/rho/cp
#
#def ierfc(z):
#    return z*erfc(z) - np.exp(-z**2)/np.sqrt(np.pi)
#
#def T_hf(F,x=0.,t=3.,k=k,kappa=kappa):
#    T = F/k * np.sqrt(4*kappa*t) * (-ierfc(x/np.sqrt(4*kappa*t)))
#    return T
#    
#F = np.linspace(0,35e6,100)
#plt.subplots()
#plt.plot(F/1e6,T_hf(F)+80, linewidth=3, label='Tile insert surface')
#plt.plot([0,200],[2623,2623], label='Mo melting point')
#plt.plot([0,200],[3422,3422], label='W melting point')
#plt.xlabel(r'q$_{\perp}$ [MW/m$^2$]', fontsize=24)
#plt.ylabel('Temperature [C]', fontsize=24)
#plt.title('Surface temperature after 3 seconds\nof constant heat flux',
#          fontsize=20)
#plt.xlim([0,np.max(F)/1e6])
##plt.ylim([1.5,11.5])
#plt.tick_params(axis='x', labelsize=18)
#plt.tick_params(axis='y', labelsize=18)
#plt.legend(loc='lower right', fontsize=18)
#plt.subplots_adjust(bottom=.14, left=.14)
#plt.show()