# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 14:10:23 2018

@author: bartonjo

Functions for file reading and writing that can be used by many 
different data analysis scripts in the LP directories.
"""

#import sys
#s = '/Users/bartonjo/PyFiles/'
#if s not in sys.path:
#    sys.path.insert(0, s)

#import os
## for stuff like:
## for file in os.listdir(path):
##        if file.endswith('.csv'):

import pandas as pd

def get_h5(filename, keys=False):
    '''
    Read in hdf5 data from storage and return it as a list of 
    dataframes.  This assumes the .h5 file has at least one pandas 
    dataframe.
    
    filename: this is the path and file name of the file to read in.
    keys: if true, the corresponding list of keys for each dataframe
    will be returned as a separate list from the data.
    
    dfs: this is the data frame list to be returned
    ks: this is the list of keys (strings) corresponding to each 
    dataframe.
    '''
    store = pd.HDFStore(filename)
    dfs,ks = [],[]
    for key in store.keys():
        print(key)
        dfs.append(pd.read_hdf(filename,key=key))
        if keys == True:
            ks.append(key[1:])  # remove the slash character from key
        
    if keys == True:
        x = [dfs, ks]
    else:
        x = dfs
    
    return x
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
