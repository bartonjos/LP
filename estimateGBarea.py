# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:12:21 2017

@author: bartonjo

Estimate the area of GB for UFG and ITER grade W
"""

import numpy as np

# Define sizes of things
rtot = 3e-3         # [m]
riter = 10e-6
rufg = 1e-6
gbthick = 1e-9      # grain boundary thickness in [m]

# Estimate GB areas
Atot = np.pi*rtot**2
Aiter_grain = np.pi*riter**2
Aufg_grain = np.pi*rufg**2
num_iter = Atot/Aiter_grain
num_ufg = Atot/Aufg_grain
def Agb(rin,rout):
    y = np.pi*(rout**2-rin**2)
    return y
Aiter_gb = Agb(riter,riter+gbthick)*num_iter
Aufg_gb = Agb(rufg,rufg+gbthick)*num_ufg

# How much more GB area does UFG W have?
R = Aufg_gb/Aiter_gb
print(R)

# Thus the amount of defects that will be created by ion damage in the
# grain boundaries will be R times more likely in UFG W.