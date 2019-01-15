# -*- coding: utf-8 -*-
"""
Script to to read and plot monitoring parameters from a tdms file produced by 
the IRSENS Labview Software.
"""

import matplotlib.pyplot as plt
import numpy as np
import tdmsProcFuncs
import copy
from scipy.stats import linregress

"""
input Data
"""

"Define file path with fiting results" 
tdms_file_path_new = 'ShortPulses/Alpha250/WorkAir#2_070918/ThroughFlow_WorkAir#2(1)_CO2-4_070918.tdms'

"""Select paremeter names to plot. Second value indicates whether the data has
to be transformed from Volt to °C or not."""

parlist = [[["p1-Cell"],[False]]]#0
parlist.append([["T1-Cell"],[False]])#1
parlist.append([["T Laser 4,3 DAQVal"],[True]])#2
parlist.append([["T Laser 7,8 DAQVal"],[True]])#3
parlist.append([["T Box DAQVal"],[True]])#4
parlist.append([["Vici Valve Pos DAQVal"],[False]])#5
parlist.append([["p Turbo DAQVal"],[False]])#6
parlist.append([["Cell Entry Valve DAQVal"],[False]])#7
parlist.append([["ByPass Valve DAQVal"],[False]])#8
parlist.append([["Turbo Valve DAQVal"],[False]])#9

"Linear detrending of parameter?"
detrend = True
detrend_par = 1 #which parameter from above to be detrended?
detrend_start = 1500 #staring line for detrending
detrend_end = -5000 #ending line for detrending

"""
Main Program
"""

"""load file if not yet loaded and make a safety copy."""
if 'tdms_file_path' not in globals() or tdms_file_path_new != tdms_file_path:
    tdms_file_path = tdms_file_path_new
    tdms_file = tdmsProcFuncs.readtdmsfile(tdms_file_path)
    tdms_file_mon = copy.deepcopy(tdms_file)

"""find group and channel names (i.e. "Fit Results" and "Line 1Amplitude, respectively)
via tdms_file.objects. Then access them as follows""" 

"""Time Stamp"""
timeStamp_mon_orig = tdms_file.object("Monitoring Vals" , "t-stamp").data
timeStamp_fit_orig = tdms_file.object("Fit Results" , "t-stamp").data

"""Transform TimeStamp into seconds after start"""
timeStamp_mon = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_orig)
timeStamp_fit = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_orig)

"""Unpack parameter values and transform to temperature where necessary"""
parvalues = []
for x in range(len(parlist)):
    parvalues.append(tdms_file.object("Monitoring Vals" , *parlist[x][0]).data)
    if parlist[x][1] == [True]:
        parvalues[x] = tdmsProcFuncs.VoltsToTempConvert(parvalues[x])

"""plot monitoring data"""
fig_1, ax_0, ax_1 = tdmsProcFuncs.plot3var_2yaxis(timeStamp_mon,parvalues[5],\
                    parvalues[1],parvalues[0],parlist[5][0],parlist[1][0],parlist[0][0])
ax_0.set_xlabel('time (sec)')
ax_0.set_ylabel('[°C]')
ax_1.set_ylabel('[atm]')
ttext = plt.title('Cell Pressure and Temperature')
plt.setp(ttext, size='large', color='r', weight='bold')

fig_2, ax_2, ax_3 = tdmsProcFuncs.plot2var_2yaxis(timeStamp_mon,parvalues[2],\
                    parvalues[3],parlist[2][0],parlist[3][0])
ax_2.set_xlabel('time (sec)')
ax_2.set_ylabel('[°C]')
ax_3.set_ylabel('[°C]')
ttext = plt.title('Laser Temperatures')
plt.setp(ttext, size='large', color='r', weight='bold')

fig_3, ax_4, ax_5 = tdmsProcFuncs.plot2var_2yaxis(timeStamp_mon,parvalues[6],\
                    parvalues[4],parlist[6][0],parlist[4][0])
ax_4.set_xlabel('time (sec)')
ax_4.set_ylabel('[? e-3 mbar ?]')
ax_5.set_ylabel('[°C]')
ttext = plt.title('Driver Temperature and Turbo Pressure')
plt.setp(ttext, size='large', color='r', weight='bold')

if detrend == True:
    timeStamp_detrend = timeStamp_mon[detrend_start:detrend_end]
    parvalues_detrend = parvalues[detrend_par][detrend_start:detrend_end]
    linreg = linregress(timeStamp_detrend,parvalues_detrend)
    plt.figure()
    plt.plot(timeStamp_detrend,parvalues_detrend-linreg.intercept,label = 'non-detrended')
    parvalues_detrend = parvalues_detrend-(linreg.slope*timeStamp_detrend+linreg.intercept)
    plt.plot(timeStamp_detrend,parvalues_detrend,label = 'detrended')
    #"plt.legend(loc='best')
