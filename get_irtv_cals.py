# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 10:58:30 2018

@author: bartonjo

Get all of the IRTV calibration shots from csv files and take the 
average intensity value for each integration time, int_times values are
given by the user. Save the temperature vs intensity (for each 
integration time) dataframe to a csv file.

With the temperature vs intensity data frame (cal) use the temp_curve 
function to fit a set of data for a particular integration time. 

With the save_fit function, use the results of this fit to convert 
measured raw data to calibrated Temperature vs time data and save it 
to a file in the directory with the rest of the experimental data. The
calibrated data will have the _Tvt.csv suffix.

Example workflow (if you have just obtained the calibration data):
1. Edit the user-defined values.

2. Run this script to define the user-defined values in memory and
initialize the functions. e.g.:
   In: runfile('get_irtv_cals.py')
   or
   In: run ./get_irtv_cals.py
   
3. In: cal = get_cal_data()

4. In: temp = np.array(cal['temp_C'])
   In: inten = np.array(cal['t0.02'])
   In: lfit, pfit, logfit = temp_curve(temp, inten)

5. Let's say you want to use the power-law fit for calibration, then:
   In: ircal = save_fit(pfit, islin=False)
   
Example if the cal data is already stored:
1. Get the calibration data from the file:
   In: cal  = pd.read_csv('SiC_cals.csv')

2. In: temp = np.array(cal['temp_C'])
   In: inten = np.array(cal['t0.02'])
   In: lfit, pfit, logfit = temp_curve(temp, inten)

3. Let's say you want to use the power-law fit for calibration, then:
   In: ircal = save_fit(pfit, islin=False, name='174508_SiC1_mid')
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
#from scipy.stats import chisquare
from scipy.stats import chi2
import matplotlib.pyplot as plt

## User defined values, if you prefer to not define them in get_cal_data
#path = '/Users/bartonjo/PyFiles/LP/irtv_images/ufg_retention_exp/cals/'
#bname = 'C_cal'  # button name used in calibration
#int_times = np.array([.2, .07, .02, .01])  # descending order!
#
## User defined values, if you prefer to not define them in save_fit
#data_path = path[:-5]  # path for the _ivt.csv data files
#data_name = '171436_U-1'  # experiment data for Tvt conversion

def get_cal_data(path, bname, int_times):
    '''
    Loop through all the files with bname (for all temperatures) 
    and get the average value of the intensity for each integration
    time.  Assume that the file name as the structure 
    temp_bname_ivt.csv (e.g. 300_W04_ivt.csv) and has two columns that
    are time_s and inten_au, where time_s is the frame * frame rate
    with the frame rate = 1 [sec/image].  Save this data frame to a
    csv file and also return it.
    '''    
    # initialize arrays
    temps = np.array([])
    itimes = np.array([])
    ierrs = np.array([])
    # create a string array for cal data dataframe column names
    dfnames = ['temp_C']
    for i in range(len(int_times)):
        dfnames.append('t' + str(int_times[i]))
        dfnames.append('e' + str(int_times[i]))
    for file in os.listdir(path):
        if file.endswith(bname + '_ivt.csv'):
            # get the temperature value from the file name
            temp = float(file[:-(len(bname + '_ivt.csv') + 1)])
            temps = np.append(temps, temp)
            # put the file in a dataframe
            df = pd.read_csv(path+file)
            # separate out intensity plateaus
            pstart, pend = separate_plateaus(df)
            # get the first plateaus and return the next row of 
            # integration times for this temperature
            itimes, ierrs = get_int_times(df,pstart,pend,
                                          int_times,itimes,ierrs)
    # create dataframe
    d = {dfnames[0]:temps}
    j=0
    for i in range(1,np.size(dfnames),2):
        d[dfnames[i]] = itimes[:,j]
        d[dfnames[i+1]] = ierrs[:,j]
        j+=1
    cal = pd.DataFrame(d)
    cal = cal.sort_values('temp_C').reset_index(drop=True)
    # save dataframe
    cal.to_csv(path+bname+'_cals.csv', index=False)
    return cal

def separate_plateaus(df):
    '''
    Separate the calibration intensities for each integration 
    time and return the start and end of each plateau.
    '''
    di = np.gradient(df.inten_au)
    pstart = np.array([0])  # plateau start indices
    pend = np.array([])  # plateau end indicies
    # to understand this loop: plot di and df.inten_au with points
    for i in range(len(di)-1):
        # find the plateau edges (they have large derivatives)
        if np.abs(di[i]) > 50 and np.abs(di[i+1]) > 50:
            if i==0:
                pstart = np.array([1])
            elif (i+1) == (len(di)-1):
                pend = np.append(pend, i)
            else:
                pend = np.append(pend, i)
                pstart = np.append(pstart, i+1)
        elif (i+1) == (len(di)-1):
            pend = np.append(pend, i+1)
    pstart = pstart.astype(int)
    pend = pend.astype(int)
    return pstart, pend                

def get_int_times(df,pstart,pend,int_times,itimes,ierrs):
    '''
    Take the longest plateaus and average over them to get the
    average intensity for that integration time. Then sort descending
    and add these intensities to the itimes matrix.
    '''
    itime = np.array([])  # integration times for a single temperature
    ierr = np.array([])  # std of measurement
    for i in range(len(int_times)):
        # take only one average for each int_times
        if (pend[0]-pstart[0]) < 10:
            avg = df.inten_au[pstart[i+1]:(pend[i+1]+1)].mean()
            er = df.inten_au[pstart[i+1]:(pend[i+1]+1)].std()
        else:
            avg = df.inten_au[pstart[i]:(pend[i]+1)].mean()
            er = df.inten_au[pstart[i]:(pend[i]+1)].std()
        itime = np.append(itime, avg)
        ierr = np.append(ierr, er)
        
    # sort itime and ierr by itime values, descending
    z = np.array([itime,ierr])
    dum = pd.DataFrame(z.transpose())
    dum = dum.sort_values(0).reset_index(drop=True)  # ascending
    itime = np.flipud(np.array(dum[0]))
    ierr = np.flipud(np.array(dum[1]))
#    print(itime)
#    print(ierr)
#    itime = np.flipud(np.sort(itime))  # sort descending
    # extend the array to add the integration times for this temp
    if np.size(itimes) == 0:
        itimes = np.append(itimes, itime)
        ierrs = np.append(ierrs, ierr)
    elif np.size(itimes) == len(int_times):
        itimes = np.append([itimes], [itime], axis=0)
        ierrs = np.append([ierrs], [ierr], axis=0)
    else:
        itimes = np.append(itimes, [itime], axis=0)
        ierrs = np.append(ierrs, [ierr], axis=0)
    return itimes, ierrs

def get_chisq(y_obs, f_exp=None, nparams=1, ers=None):
    '''
    Calculate the reduced chi squared with uniform errors.

    y_obs: data for fitting.
    f_exp: value of the fit corresponding to y_obs. Default will be the
        mean of y_obs if no fit data is given.
    nparams: number of parameters in fitting equation.

    Returns xsq and pval
    xsq: reduced chi square value.
    pval: p-value for the distribution.    
    '''
    if f_exp is None:
        f_exp = np.atleast_1d(y_obs.mean())
    else:
        # remove f_exp = 0 values
        y_obs = np.array([y_obs[i] for i in range(np.size(f_exp)) \
                if f_exp[i] != 0])
        f_exp = np.array([f_exp[i] for i in range(np.size(f_exp)) \
                if f_exp[i] != 0])
    
    # degrees of freedom
    dof = np.size(y_obs)-nparams
    
    # Pearson chi-squared calculation:
    terms = (y_obs - f_exp)**2/f_exp
    
    # reduced chi squared
    xsq = terms.sum()/dof
    # p-value is the sf = 1 - cdf of the chi squared distribution
    pval = chi2.sf(xsq,dof)
    
    return xsq, pval
    
    
def temp_curve(temp, inten, ers=None):
    '''
    Input numpy arrays of temperature and intensities and perform
    a linear and power-law fit. Return the coefficients of the 
    fitting functions.
    The intensity array you supply is for just one integration time!
    Choose the intensities with integration time used for the 
    experimental data.
    '''
    
    # perform a linear least squares fit
    z = np.polyfit(inten, temp, 1, full=True)
    res = float(z[1])  # residuals of fit: sum((y-f)**2)
    rsq = 1-res/np.sum((temp-temp.mean())**2)  # R^2 value of fit
    linfit = z[0]  # fit coefficients
    p = np.poly1d(z[0])  # linear function
    xsl, pl = get_chisq(temp, f_exp=p(inten), nparams=2)
    print('Reduced chi-sq value for linear fit = %2.4f' % xsl)

    # perform a fit with a power-law function
    try:
        si = inten-np.min(inten)
        st = temp-np.min(temp)
        popt, pcov = curve_fit(plaw2, si, st)
        xsp, pp = get_chisq(st, f_exp=plaw2(si,*popt),
                            nparams=2)
        popt = np.array([popt[0],-np.min(inten),popt[1],np.min(temp)])
#            popt, pcov = curve_fit(plaw, inten, temp)
#            xsp, pp = get_chisq(temp,f_exp=plaw(inten,*popt),nparams=4)
        print('Reduced chi-sq value for power-law fit = %2.4f' % xsp)
    except:
        print('Power-law fit could not converge. Try another fit.')
        popt = 0
        pass

    # perform a fit of intensity vs T, then invert
    try:
        popt2, pcov2 = curve_fit(rev, temp, inten)    
        xsr, pr = get_chisq(temp, f_exp=rev_fit(inten,*popt2), 
                            nparams=3)
        print('Reduced chi-sq value for log fit = %2.4f' % xsr)
    except:
        print('Log fit could not converge. Try another fit.')
        popt2 = 0
        pass

    # plotting
    plt.figure()    
    plt.plot(inten,temp,'s',label='measured')
    if ers is not None:
        plt.errorbar(inten, temp, yerr=ers, fmt=' b')
    # visualize linear fit
    lstr = r'T = %5.2f $\times$ I + %5.2f' % tuple(linfit) 
    x = np.linspace(np.min(inten)*.75, np.max(inten)*1.25,200)
    plt.plot(x,p(x),label=lstr)
    erst = 'Linear fit: \n' r'R$^2$ = %0.4f' '\n' % rsq
    erst = erst + r'$\chi^2_{\nu}$ = %2.4f' %xsl + \
            '\n' 'p-value = %0.4f' % pl
    plt.text(np.min(x),np.max(temp), erst)
    # visualize power-law fit
    if 'pp' in locals():
        pstr = 'T = ' + \
          r'%5.2f $\times$ (I + %7.2f)$^{%4.2f}$ + %5.2f' % tuple(popt)        
        plt.plot(x,plaw(x,*popt),label=pstr) 
        erstp = 'Power-law fit: \n' r'$\chi^2_{\nu}$ = %2.4f' % xsp + \
                '\n p-value = %0.4f' % pp
        plt.text(np.min(x), np.max(temp)*.5, erstp)
    # visualize reverse power (log) fit
    if 'pr' in locals():
        a,b,c = popt2
        rstr = r'T = $\log_{%0.4f}$(I-%0.1f) - %0.1f' % tuple([a,c,b])
        plt.plot(x,rev_fit(x,*popt2), label=rstr)
        erstr = 'Log fit: \n' r'$\chi^2_{\nu}$ = %2.4f' % xsr + \
                    '\n p-value = %0.4f' % pr
        plt.text(np.min(x)*1.5, np.max(temp), erstr)
    # edit plot
    plt.tick_params(axis='x', labelsize=18)
    plt.tick_params(axis='y', labelsize=18)
    plt.ylabel('Temperature [C]', fontsize=22)
    plt.xlabel('Intensity [au]', fontsize=22)
    plt.legend(loc='lower right', fontsize=12)
    plt.subplots_adjust(bottom=.11, left=.14)
    plt.show()
    return linfit, popt, popt2

def plaw2(x,a,c):
    '''
    power law function relating T vs intensity. 
    note, it is assumed 0 < c <= 1 and (x,y) starts at the origin
    '''
    return a * (x)**c
    
def plaw(x,a,b,c,d):
    '''
    power law function relating T vs intensity. 
    note, it is assumed 0 < c <= 1
    '''
    return a * (x + b)**c + d
    
def rev(x, a, b, c):
    '''
    fuction to fit intensity vs T
    '''
    return a**(x + b) + c

def rev_fit(x, a, b, c):
    '''
    function of the inversion of rev to give T vs intensity
    '''
    return np.log(x-c)/np.log(a) - b

def save_fit(coeffs, path, name, islin=True,logfit=False):
    '''
    Specify which fit you want to use and use that calibration on a
    given file.  Save the temperature calibrated file and also return
    it as a dataframe.
    '''
    # import data from _ivt.csv file into a dataframe
    fname = name + '_ivt.csv'
    uncal = pd.read_csv(path+fname)
    
    if islin == True:
        # get linear function from coeffs 
        p = np.poly1d(coeffs)
        # convert intensity to temperature via the p() function
        temp = p(np.array(uncal.inten_au))
    elif logfit == False:
        # convert intensity to temperature via the plaw() function
        temp = plaw(np.array(uncal.inten_au),*coeffs)
    else:
        # convert intensity to temperature via the rev_fit() function
        temp = rev_fit(np.array(uncal.inten_au),*coeffs)
        
    # create a dataframe with time_s and Temp_C columns
    ircal = pd.DataFrame({'time_s':np.array(uncal.time_s),
                          'temp_C':temp})
    # save the dataframe in the data_path directory
    calfname = name + '_Tvt.csv'
    ircal.to_csv(path + calfname, index=False)
    
    return ircal





















    