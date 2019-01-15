# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 14:40:58 2017

File that contains all packages that need to be importad and all the 
definitions needed for processing the data of the LabView IRSENS Software
stored in tdms file format.

@author: beb
"""

"""
Imports
"""

import nptdms
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import copy
import allantools
from astropy.stats import sigma_clip

"""
functions
"""
def readtdmsfile(filepath):
    """ read tdms file at filepath """
    tdms_file = nptdms.TdmsFile(filepath)
    return tdms_file

def VoltsToTempConvert(voltsdata):
    """converts volts measured at a temperature sensor into temperature"""
    temp = 1/(1.1279*10**(-3)+2.3429*10**(-4)*np.log(voltsdata/0.0001)\
              +8.7298*10**(-8)*(np.log(voltsdata/0.0001))**3)-273.15
    return temp

def plot1var(x,y,label):
    """plots x vs. y"""
    fig, ax = plt.subplots()
    plt1 = ax.plot(x,y,'r',label=label)
    ax.legend(loc='best')
    ax.grid(True)
    return fig, ax, plt1

def plot2var_2yaxis(x,y1,y2,label1,label2):
    """plots x vs. y1 and y2 using two independent y-axis on each side of the plot"""
    fig, ax_1 = plt.subplots()
    ax_1.plot(x,y1,'g',label=label1)
    ax_1.legend(loc='lower left')
    ax_1.grid(True)
    ax_2 = ax_1.twinx()
    ax_2.plot(x,y2,'m',label=label2)
    ax_2.legend(loc='lower right')
    return fig, ax_1, ax_2
    
def plot2var_1yaxis(x,y1,y2,label1,label2):
    """plots x vs. y1 and y2 using the same y-axis"""
    fig, ax = plt.subplots()
    ax.plot(x,y1,'c',label=label1)
    ax.plot(x,y2,'y',label=label2)
    ax.legend(loc='best')
    ax.grid(True)
    return fig, ax

def plot2var_1yaxis_xaxisDiff(x1,y1,x2,y2,label1,label2):
    """plots x1 vs. y1 and x2 vs. y2 using the same y-axis."""
    plt.figure()
    plt.hold(True)
    ax_1 = plt.plot(x1,y1,'r',label=label1)
    ax_1 = plt.plot(x2,y2,'b',label=label2)
    plt.grid(True)
    plt.legend(loc='best')
    plt.hold(False)
    plt.show()
    return ax_1

def plot3var_2yaxis(x,y1,y2,y3,label1,label2,label3):
    """plots x vs. y1 and y2 on the same y-axis and x vs. y3 on a seperate y-axis"""
    fig, ax_1 = plt.subplots()
    ax_1.plot(x,y1,'r',label=label1)
    ax_1.plot(x,y2,'b',label=label2)
    ax_1.legend(loc='lower left')
    ax_1.grid(True)
    ax_2 = ax_1.twinx()
    ax_2.plot(x,y3,'k',label=label3)
    ax_2.legend(loc='lower right')
    return fig, ax_1, ax_2

def plot4var_2yaxis(x,y1,y2,y3,y4,label1,label2,label3,label4):
    """plots x vs. y1 and y2 on the same y-axis and x vs. y3 and y4 on a seperate y-axis"""
    fig, ax_1 = plt.subplots()
    ax_1.plot(x,y1,'r',label=label1)
    ax_1.plot(x,y2,'b',label=label2)
    ax_1.legend(loc='lower left')
    ax_1.grid(True)
    ax_2 = ax_1.twinx()
    ax_2.plot(x,y3,'k',label=label3)
    ax_2.plot(x,y4,'m',label=label4)
    ax_2.legend(loc='lower right')
    return fig, ax_1, ax_2
    
def plotallan_adev(plt, y, rate, taus, style):
    (t2, ad, ade, adn) = allantools.adev(y, rate=rate, data_type="freq", taus=taus)
    plt.loglog(t2, ad, style)
    return t2, ad, ade, adn

def plotallan_adev_phase(plt, y, rate, taus, style):
    (t2, ad, ade, adn) = allantools.adev(y, rate=rate, taus=taus)
    plt.loglog(t2, ad, style)
    return t2, ad, ade, adn

def plotallan_oadev(plt, y, rate, taus, style):
    (t2, ad, ade, adn) = allantools.oadev(y, rate=rate, data_type="freq", taus=taus)
    plt.loglog(t2, ad, style)
    ad_min_y = round(ad.min(),4)#save Allan minimum
    ad_min_x = round(t2[np.argmin(ad)])#save corresponding averaging time/tau
    string = "Allan Minimum " + str(ad_min_y) + "\nIntegration time " + str(ad_min_x)
    plt.text(ad_min_x, ad_min_y+(ad.max()-ad_min_y)/2, string)
    plt.ylabel('Allan Deviation [per mil]')
    plt.xlabel('Time [sec]')
    ttext = plt.title('Allan Deviation')
    plt.setp(ttext, size='large', color='r', weight='bold')
    return t2, ad, ade, adn

def plotallan_oadev_phase(plt, y, rate, taus, style):
    (t2, ad, ade, adn) = allantools.oadev(y, rate=rate, taus=taus)
    plt.loglog(t2, ad, style)
    return t2, ad, ade, adn
   
def plotallan_mdev(plt, y, rate, taus, style, label=""):
    (t2, ad, ade, adn) = allantools.mdev(y, data_type='freq', rate=rate, taus=taus)
    plt.loglog(t2, ad, style,label=label)
    return t2, ad, ade, adn

def plotallan_mdev_phase(plt, y, rate, taus, style, label="",alpha=1.0):
    (t2, ad, ade, adn) = allantools.mdev(y,data_type='phase', rate=rate, taus=taus)
    plt.loglog(t2, ad, style, label=label,alpha=alpha)
    return t2, ad, ade, adn
    
def plotline(plt, alpha, start, taus, style):
    """ plot a line with the slope alpha """
    y = [start*pow(tt, alpha) for tt in taus]
    plt.loglog(taus, y, style)
    
def TimeStampTransform(timeStamp):
    """Transform timeStamp into seconds after start"""
    timeStamp = mdates.date2num(timeStamp)
    timeStamp = timeStamp-timeStamp[1] #This is time expressed in days
    timeStamp = timeStamp*3600*24 #Convert to Seconds 
    return timeStamp

def outlierFilt(data,sigma):
    """Find outliers and replace them with subsequent values"""
    datacopy = copy.copy(data)
    flags = sigma_clip(datacopy, sigma=sigma).mask
    flagratio = flags.sum()/len(datacopy)
    flags[-1] = False #Prevent error in case the last element is flagged
    flags_shift = np.roll(flags,1) #shift the flags by one element
    for x in range(5): #loop over replacement in case neighbours are flagged
        datacopy[flags] = datacopy[flags_shift]
    return datacopy, flagratio

"""fast way to find nearest element from left in an array - Josh Albert, Stack Overflow"""
def bisection(array,value):
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl
    