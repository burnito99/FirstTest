# -*- coding: utf-8 -*-
"""
Script to to read CH4 and N2O fitting data and calculate the basic parameters 
based on the tdms file produced by the IRSENS Labview Software.
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import copy
import tdmsProcFuncs

"""
input Data
"""
"Define file path with fiting results" 
tdms_file_path_new = 'ShortPulses/Redpitaya/StandAir_080617/StandAir_CH4_080617.tdms'

"Foreground partial pressure fit parameter of N2O and CH4 line"
line_N2O = "Line 0Amplitude"
line_CH4 = "Line 1Amplitude"
label_CH4 = 'CH4'
label_N2O = 'N2O'

"N2O concentration of N2O set in Hitran [ppb]"
pN2O = 100
"Cell pressure set in Hitran [atm]"
cellPress = 0.005
"CH4 concentration of CH4 set in Hitran [ppb]"
pCH4 = 1000

"other important parameters saved within 'Fit Results' (e.g. baseline parameters) "
parlist_ch4 = [["P-Coeff 0"]]
parlist_ch4.append(["P-Coeff 1"])
parlist_ch4.append(["P-Coeff 2"])
parlist_ch4.append(["P-Coeff 3"])
parlist_ch4.append(["Zero Level"])

"""
Main Program
"""

"""load file if not yet loaded and make a safety copy."""
if 'tdms_file_path' not in globals() or tdms_file_path_new != tdms_file_path:
    tdms_file_path = tdms_file_path_new
    tdms_file = tdmsProcFuncs.readtdmsfile(tdms_file_path)
    tdms_file_CH4 = copy.deepcopy(tdms_file)

"""find group and channel names (i.e. "Fit Results" and "Line 1Amplitude, respectively)
via tdms_file.objects. Then access them as follows""" 

"""Time Stamp"""
timeStamp_fit_orig = tdms_file.object("Fit Results" , "t-stamp").data
timeStamp_mon_orig = tdms_file.object("Monitoring Vals" , "t-stamp").data

"""Transform TimeStamp into seconds after start"""
timeStamp_fit = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_orig)
timeStamp_mon = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_orig)

"""Check if time stamp of Monitoring Vals and Fit Results fit"""
if len(timeStamp_fit) != len(timeStamp_mon):
    warnings.warn('Time stamp of monitoring and fit values not same length!')
else:
    if not np.allclose(timeStamp_fit,timeStamp_mon, atol = 1, rtol = 0):
        warnings.warn('Time stamp of monitoring and fit values is shifted > 1sec.') 

"""Amplitudes of main lines"""
lineAmp_N2O = tdms_file.object("Fit Results" , line_N2O).data
lineAmp_CH4 = tdms_file.object("Fit Results" , line_CH4).data

"""Cell Pressure"""
p1Cell = tdms_file.object("Monitoring Vals" , "p1-Cell").data
p1Cell = np.mean(p1Cell) #remove "noisy" pressure reading

"""Calculate concentrations of N2O and CH4"""
lineCons_N2O = lineAmp_N2O*pN2O*cellPress/p1Cell
lineCons_CH4 = lineAmp_CH4*pCH4*cellPress/p1Cell

"""Unpack other parameters"""
parvalues_ch4 = []
for x in range(len(parlist_ch4)):
    parvalues_ch4.append(tdms_file.object("Fit Results" , *parlist_ch4[x]).data)

"""plot N2O and CH4 data"""
fig_ch4, ax_ch4, ax_n2o = tdmsProcFuncs.plot2var_2yaxis(timeStamp_fit,lineCons_CH4,\
                    lineCons_N2O,label_CH4,label_N2O)
ax_ch4.set_ylabel('CH4 [ppb]')
ax_ch4.set_xlabel('time (sec)')
ax_n2o.set_ylabel('N2O [ppb]')
ttext = plt.title('CH4 and N2O in cell')
plt.setp(ttext, size='large', color='r', weight='bold')

"""plot other data"""
fig_pol, ax_0, ax_1 = tdmsProcFuncs.plot4var_2yaxis(timeStamp_fit,parvalues_ch4[0],\
                    parvalues_ch4[1],parvalues_ch4[2],parvalues_ch4[3],parlist_ch4[0],\
                    parlist_ch4[1],parlist_ch4[2],parlist_ch4[3])
ax_0.set_xlabel('time (sec)')
ax_0.set_ylabel('y1 and y2 values')
ax_1.set_ylabel('y3 and y4 values')
ttext = plt.title('Baseline Polynom Parameters of 7.8 um laser')
plt.setp(ttext, size='large', color='r', weight='bold')
