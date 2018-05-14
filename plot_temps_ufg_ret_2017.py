# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:21:54 2018

@author: bartonjo

Scripts to plot IRTV and TC data for the retention experiment.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
s = '/Users/bartonjo/PyFiles/'
if s not in sys.path:
    sys.path.insert(0, s)
s = '/Users/bartonjo/PyFiles/LP/'
if s not in sys.path:
    sys.path.insert(0, s)
from get_data import get_h5

path = '/Users/bartonjo/PyFiles/LP/irtv_images/ufg_retention_exp/'

buttons = {0:'I-7', 1:'I-3', 2:'I-4', 
           3:'U-5', 4:'U-1', 5:'U-2', 6:'U-4',
           7:'C_i', 8:'C_m', 9:'C_o', 10:'C_m2'}

props = {'I-7':'ITER W, 0.6 dpa, inner', 
         'I-3':'ITER W, 0.6 dpa, middle',
         'I-4':'ITER W, 0 dpa, outer', 'U-5':'UFG W, 0.6 dpa, inner', 
         'U-1':'UFG W, 0 dpa, middle', 'U-2':'UFG W, > 0.6 dpa, outer', 
         'U-4':'UFG W, 0.06 dpa, middle', 'C_i':'cap downstream of E', 
         'C_m':'cap downstream of G', 'C_o':'cap downstream of F', 
         'C_m2':'cap downstream of D'}

def data_list(bname, tc_temp=False):
    '''
    return a list of dataframes for all shots for a button.
    
    arguments
    bname: string value that is the button's name, like 'U-1'
    tc_temp: if true, then get the TC temperature data for 'U-1'
    
    returns
    dfs: dataframe list with each element corresponding to subsequent 
    plasma discharges
    '''
    
    # switch i-4 and i-3 temperature data
    if bname == 'I-3':
        bname = 'I-4'
    elif bname == 'I-4':
        bname = 'I-3'
        
    dfs = []
    if tc_temp == True:
        dfs = get_h5(path + 'tc_temps/tc_data.h5')
    else:
        fname = '171436_' + bname + '_Tvt.csv'
        for i in range(10):
            dfs.append(pd.read_csv(path+fname))
            j = int(fname[4:6])+1
            fname = fname[:4] + str(j) + fname[6:]
    
    return dfs

## ---------------------------------------------------------------------
## calibrate all buttons
## ---------------------------------------------------------------------
#
#import get_irtv_cals as gc
#
#cal  = pd.read_csv(path + 'cals/W_cal_cals.csv')
##cal  = pd.read_csv(path + 'cals/C_cal_cals.csv')
#
## all buttons of first shot have different intgration time
#temp = np.array(cal['temp_C'])
#inten = np.array(cal['t0.02'])
#lfit, pfit, logfit = gc.temp_curve(temp, inten)
#for i in range(7):
#    name='171436_'+buttons[i]
#    ircal = gc.save_fit(logfit, path, name, islin=False, logfit=True)
##for i in [10]: #range(7,10):
##    name='171436_'+buttons[i]
##    ircal = gc.save_fit(logfit, path, name, islin=False, logfit=True)
#
## all buttons for the remaining 9 shots have same integration time
#inten = np.array(cal['t0.2'])
#lfit, pfit, logfit = gc.temp_curve(temp, inten)
#for i in range(7):
#    fname = '171437_' + buttons[i]
##for i in [10]: # range(7,10):
##    fname = '171437_' + buttons[i]
#    for j in range(1,10):
#        print(fname)
#        ircal = gc.save_fit(logfit, path, fname, 
#                            islin=False, logfit=True)
#        k = int(fname[4:6])+1
#        fname = fname[:4] + str(k) + fname[6:]

## ---------------------------------------------------------------------
## Plot IRTV data for shot 171444 for UFG and ITER buttons
## ---------------------------------------------------------------------
#i7 = pd.read_csv(path + '171444_I-7_Tvt.csv')
#i3 = pd.read_csv(path + '171444_I-3_Tvt.csv')
#i4 = pd.read_csv(path + '171444_I-4_Tvt.csv')
#
#plt.figure()
#
#plt.plot(i4.time_s,i4.temp_C,label='W middle')
#plt.plot(i3.time_s,i3.temp_C,label='W outer')
#plt.plot(i7.time_s,i7.temp_C,label='W inner')
#plt.xlim([-.5,6])
#plt.ylim([60,300])
#plt.legend(loc='best')
#plt.title('171444 ITER W buttons ir data')
#
#u5 = pd.read_csv(path + '171444_U-5_Tvt.csv')
#u1 = pd.read_csv(path + '171444_U-1_Tvt.csv')
#u2 = pd.read_csv(path + '171444_U-2_Tvt.csv')
#u4 = pd.read_csv(path + '171444_U-4_Tvt.csv')
#
#plt.figure()
#
#plt.plot(u1.time_s,u1.temp_C,label='U middle')
##plt.plot(u4.time_s,u4.temp_C,label='U middle 2')
#plt.plot(u2.time_s,u2.temp_C,label='U outer')
#plt.plot(u5.time_s,u5.temp_C,label='U inner')
#plt.xlim([-.5,6])
#plt.ylim([60,300])
#plt.legend(loc='best')
#plt.title('171444 UFG W buttons ir data')



## ---------------------------------------------------------------------
## Plot all shots for a particular button
## ---------------------------------------------------------------------
#
#bname = 'C_o'
#dfs = data_list(bname)
##tcs = data_list('U-1',tc_temp=True)
#plt.figure()
#for i in range(10):
#    plt.plot(dfs[i].time_s,dfs[i].temp_C,label=str(171436 + i))
##    plt.plot(tcs[i].time/1e3,tcs[i].temp,label=str(171436 + i))
#
#plt.xlim([-.5,6])
#plt.ylim([60,300])
#plt.legend(loc='best')
#plt.title('All IR temperature data for button ' + bname)





# ---------------------------------------------------------------------
# Plot peak temperatures at end of shot vs shot # to show ratcheting
# ---------------------------------------------------------------------

## single button
#bname = 'U-2'
#dfs = data_list(bname)
#
#m = []
#for i in range(10):
#    z = np.min(np.where(dfs[i].time_s > 5.1))
#    m.append(dfs[i].temp_C[:z].max())
#plt.figure()
#plt.plot(np.array(range(10))+1, m, '-s', label=bname)

# multiple buttons
plt.figure()
fsize=22  # font size of axis text
sym = {0:'o', 1:'^', 2:'D', 3:'s', 4:'x', 5:'*', 6:'>'}
for j in range(7):
    dfs = data_list(buttons[j])
    m = []
    for i in range(10):
        z = np.min(np.where(dfs[i].time_s > 5.1))
        m.append(dfs[i].temp_C[:z].max())
    plt.plot(np.array(range(10)), m, '-'+sym[j], markersize=13, 
             label=props[buttons[j]])

plt.xlim([-.5,9.5])
plt.ylim([100,280])
plt.xlabel('Shot # from 171436', fontsize=fsize)
plt.ylabel('Max Temperature [C]', fontsize=fsize)
plt.tick_params(axis='both', labelsize=fsize-2)
plt.legend(bbox_to_anchor=(.5,1.2), fontsize=12)
plt.subplots_adjust(top=.85, bottom=.12, left=.13)
plt.show()


## ---------------------------------------------------------------------
## Plot temperature histories for a single shot, comparing buttons
## ---------------------------------------------------------------------
#
#shots = 171437 + np.array(range(9))
#
#plt.figure()
#fsize=22  # font size of axis text
#for i in [2,4]:  # do this for every button
#    dfs = data_list(buttons[i])
#    for j in (shots-171436):  # do this for every shot you want
##        if np.size(shots) > 1:
##            lab = props[buttons[i]] + ', ' + str(j + 171436)
##        else:
##            lab = props[buttons[i]]
##        plt.plot(dfs[j].time_s, dfs[j].temp_C, label=lab)
#        if i == 2:
#            plt.plot(dfs[j].time_s, dfs[j].temp_C, '-r')
#        else:
#            plt.plot(dfs[j].time_s, dfs[j].temp_C, '-k')
#
#if np.size(shots) == 1:
#    plt.title(str(shots[0]))
#plt.plot([1.1,1.1],[40,220],'k')
#plt.plot([5.1,5.1],[40,220],'k')
#plt.xlabel('t [s]', fontsize=fsize)
#plt.ylabel('Temperature [C]', fontsize=fsize)
#plt.tick_params(labelsize=fsize-2)
##plt.legend(loc='best', fontsize=12)
#plt.subplots_adjust(bottom=.12)
#plt.show()

















