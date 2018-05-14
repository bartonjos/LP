# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:21:54 2018

@author: bartonjo

Scripts to plot IRTV and TC data for the heat flux experiment.
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

path = '/Users/bartonjo/PyFiles/LP/irtv_images/ufg_heatflux_exp/'

buttons = {0:'SiC5', 1:'cap', 
           2:'SiC1', 3:'SiC1_low', 4:'SiC1_mid', 5:'SiC1_top', 
           6:'SiC2', 7:'SiC2_low', 8:'SiC2_mid', 9:'SiC2_top',
           10:'SiC6',
           11:'TSC1', 12:'TSC1_low', 13:'TSC1_mid', 14:'TSC1_top',
           15:'TSC2_high', 16:'TSC2_low', 17:'TSC2', 18:'TSC2_top',
           19:'TSC4', 20:'UFGNi_W2'}

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
        
    dfs = []
    if tc_temp == True:
        # to do: add temperature data 
        dfs = get_h5(path + 'tc_temps/tc_data.h5')  #!!!!! wrong
    else:
        fname = '171436_' + bname + '_Tvt.csv'
        for i in range(10):
            dfs.append(pd.read_csv(path+fname))
            j = int(fname[4:6])+1
            fname = fname[:4] + str(j) + fname[6:]
    
    return dfs

# ---------------------------------------------------------------------
# calibrate all buttons
# ---------------------------------------------------------------------

import get_irtv_cals as gc

# which calibration?
#cal  = pd.read_csv(path + 'cals/W04_cals.csv')
cal  = pd.read_csv(path + 'cals/SiC_cals.csv')
#cal  = pd.read_csv(path + 'cals/maxTi_cals.csv')

# which shot range?
#r = [x for x in range(8,16) if x!=11]
r = range(1,7)

# which buttons?
# L-mode shots:----------- 8-15
# max phase buttons 11-19
# SiC buttons 2-10
# ufg-ni button 20
# H-mode shots:----------- 1-6
# SiC5 button 0
# cap button 1
b = range(1)

# all buttons have same intgration time
temp = np.array(cal['temp_C'])
inten = np.array(cal['t0.02'])
lfit, pfit, logfit = gc.temp_curve(temp, inten)

# which fit? These are arguments in save_fit()
f = lfit
lin = True
lf = False

# Apply fit to all appropriate buttons for all shots
for i in b:  #cycle over buttons
    fname = '1745##_' + buttons[i]
    for j in r:   # cycle over shots
        if j < 10:
            fname = fname[:4] + '0' + str(j) + fname[6:]
        else:
            fname = fname[:4] + str(j) + fname[6:]
        print(fname)
        ircal = gc.save_fit(f, path, fname, islin=lin, logfit=lf)


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



# ---------------------------------------------------------------------
# Plot all shots for a particular button
# ---------------------------------------------------------------------

bname = 'C_o'
dfs = data_list(bname)
#tcs = data_list('U-1',tc_temp=True)
plt.figure()
for i in range(10):
    plt.plot(dfs[i].time_s,dfs[i].temp_C,label=str(171436 + i))
#    plt.plot(tcs[i].time/1e3,tcs[i].temp,label=str(171436 + i))

plt.xlim([-.5,6])
plt.ylim([60,300])
plt.legend(loc='best')
plt.title('All IR temperature data for button ' + bname)





## ---------------------------------------------------------------------
## Plot peak temperatures at end of shot vs shot # to show ratcheting
## ---------------------------------------------------------------------
#
### single button
##bname = 'U-2'
##dfs = data_list(bname)
##
##m = []
##for i in range(10):
##    z = np.min(np.where(dfs[i].time_s > 5.1))
##    m.append(dfs[i].temp_C[:z].max())
##plt.figure()
##plt.plot(np.array(range(10))+1, m, '-s', label=bname)
#
## multiple buttons
#plt.figure()
#for j in [7,9,10]:
#    dfs = data_list(buttons[j])
#    m = []
#    for i in range(10):
#        z = np.min(np.where(dfs[i].time_s > 5.1))
#        m.append(dfs[i].temp_C[:z].max())
#    plt.plot(np.array(range(10))+1, m, '-s', label=buttons[j])
#
#plt.xlim([.5,10.5])
#plt.ylim([100,250])
#plt.xlabel('Shot', fontsize=18)
#plt.ylabel('Max Temperature [C]', fontsize=18)
#plt.tick_params(axis='both', labelsize=18)
#plt.legend(loc='best')
#plt.show()






















