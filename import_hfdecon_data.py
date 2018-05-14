# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 16:38:45 2017

@author: bartonjo

Import data from hfdecon_... files
The assumption is that this script is run in the ~/PyFiles/LP directory
"""

import numpy as np

# data attributes from files
shot = 167408
tcname = 10
pname = 19

def getdecondata(shot=shot,tcname=tcname,pname=pname):
    # assemble file names
    ftc = 'decondata/'+'hfdecon_'+str(shot)+'_'+str(tcname)+'tc.dat'
    flp = 'decondata/'+'hfdecon_'+str(shot)+'_'+str(pname)+'p.dat'
    fir = 'decondata/'+'hfdecon_'+str(shot)+'_ir.dat'
    fe = 'decondata/'+'hfdecon_'+str(shot)+'_elmf.dat'
    fts = 'decondata/'+'hfdecon_'+str(shot)+'_dts.dat'
    ffs3 = 'decondata/'+'hfdecon_'+str(shot)+'_fs3.dat'
    ffs5 = 'decondata/'+'hfdecon_'+str(shot)+'_fs5.dat'
    flpmed = 'decondata/'+'hfdecon_'+str(shot)+'_'+\
            str(pname)+'pmed.dat'    
    fstart = 'decondata/'+'hfdecon_'+str(shot)+'_elmstart.dat'
    fosp = 'decondata/'+'hfdecon_'+str(shot)+'_osp.dat'
    fdw = 'decondata/'+'hfdecon_'+str(shot)+'_dwmhdf.dat'    
    
    # get data from file and put into np arrays
    tct, tcTemp, tcderiv = np.loadtxt(ftc, unpack=True)
    lpt, jsat, Te, ang, lpq = np.loadtxt(flp, unpack=True)
    try:
        irt, irq, irq2 = np.loadtxt(fir, unpack=True)
    except ValueError:
        print('Only (r) integrated irq data for this shot')
        irt, irq = np.loadtxt(fir, unpack=True)
    dsselmt, dsselmf = np.loadtxt(fe, unpack=True)
    tst, tsq, tsq1 = np.loadtxt(fts, unpack=True)
    try:
        fs3t, fs3 = np.loadtxt(ffs3, unpack=True)
        fs5t, fs5 = np.loadtxt(ffs5, unpack=True)
    except FileNotFoundError:
        print('No filterscope data for this shot')
        pass
    lptmed, lpqmed = np.loadtxt(flpmed, unpack=True)
    try:
        elmstart, esize, dur = np.loadtxt(fstart, unpack=True)
    except:
        print('No elm data for this shot')
        pass
    try:
        ospt, osp = np.loadtxt(fosp, unpack=True)
    except FileNotFoundError:
        print('No OSP data for this shot')
        pass
    try:
        dwt, dw = np.loadtxt(fdw, unpack=True)
    except:
        print('No dwmhd data for this shot')
        pass
    
    # assign this data to a new class object
    class shotdata(object):
        """represents data from the shot"""
        
    s = shotdata()
    s.tct = tct
    s.tcTemp = tcTemp
    s.tcderiv = tcderiv
    s.lpt = lpt
    s.jsat = jsat
    s.Te = Te
    s.ang = ang
    s.lpq = lpq
    s.irt = irt
    s.irq = irq
    try:
        s.irq2 = irq2
    except:
        pass
    s.dsselmt = dsselmt
    s.dsselmf = dsselmf
    s.tst = tst
    s.tsq = tsq
    s.tsq1 = tsq1
    try:
        s.fs3t = fs3t
        s.fs3 = fs3
        s.fs5t = fs5t
        s.fs5 = fs5 
    except:
        pass
    s.lptmed = lptmed
    s.lpqmed = lpqmed
    try:
        s.elmstart = elmstart
        s.esize = esize
        s.dur = dur
    except:
        pass
    try:
        s.ospt = ospt
        s.osp = osp
    except:
        pass
    try:
        s.dwt = dwt
        s.dw = dw
    except:
        pass
    
    return s

if __name__ == '__main__':
    # import_hfdecon_data.py executed as script
    # do something
    s = getdecondata(shot=shot,tcname=tcname,pname=pname)



















