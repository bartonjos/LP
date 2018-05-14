# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:32:01 2018

@author: bartonjo

Get ascii data from saved irtv images on DiMES and produce pixel 
intensity vs time data. Save to csv file.

For calibration files, extend the path to "/cals/" and change the 
frame rate to 1 (i.e. frate=1.).  This assumes you saved the 
calibration files into a subdirectory called cals.

Example to run multiple files:
1) comment out fname in the Inputs section.
2) run the following commands in the console:
fname = '171436_I-3_ir_image.dat'
for i in range(10):
    print(fname)
    runfile('get_irtv_imageVtime.py')
    j = int(fname[4:6])+1
    fname = fname[:4] + str(j) + fname[6:]
"""

#import sys
#s = '/Users/bartonjo/PyFiles/LP'
#if s not in sys.path:
#    sys.path.insert(0, s)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from cookbook import savitzky_golay as smooth

# Inputs unique to this dataset----------------------------------------
path = '/Users/bartonjo/PyFiles/LP/irtv_images/ufg_heatflux_exp/'
fname = '174508_SiC2_low_ir_image.dat'
frate = 1./358.1#1./126.33 ret. #1./358.1 hf.  # frame rate [s/fr]
sname = fname[:-13] + '_ivt.csv'  # csv file name to be created

# Get ascii data from file---------------------------------------------
ncols = np.size(pd.read_table(path+fname, sep='\s+', 
                              skiprows=[0]).columns)
colnames = ['frame']
for i in range(1,ncols):
    colnames.append('p'+str(i))
raw = pd.read_table(path+fname, sep='\s+', skiprows=[0], 
                    names=colnames)

# Convert frames to time-----------------------------------------------
# find how many lines in each frame and how many frames there are
nlines = raw[raw.frame == 0].shape[0]
nframes = int(float(raw.shape[0])/float(nlines))
# use the total number of lines in raw to create the time array
t = np.array(range(nframes))*frate

# Get average pixel intensity for each frame---------------------------
inten = np.array([])
for i in range(nframes):
#    inten = np.append(inten,
#                      np.max(np.max(raw.iloc[(5*i):(5+5*i),1:ncols])))
    inten = np.append(inten,np.mean(
                            np.mean(raw.iloc[(5*i):(5+5*i),1:ncols])))


# Using the time and intensity arrays, create a dataframe--------------
ivt = pd.DataFrame({'time_s':t, 'inten_au':inten})

plt.plot(ivt.time_s, ivt.inten_au, '-*')

# Save dataframe to csv file-------------------------------------------
ivt.to_csv(path+sname, index=False)












