# -*- coding: utf-8 -*-
"""
Script to to read CO2 fitting data and calculate the basic parameters based on 
the tdms file produced by the IRSENS Labview Software. Comparted to 
DataProcess_CO2.py, this script takes the results from different files that 
contain the data from different acquisition windows
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import copy
import tdmsProcFuncs
from astropy.stats import sigma_clip
import allantools

"""
input Data
"""

"Foreground partial pressure fit parameter of main CH4 line"
line_CH4 = "Line 0Amplitude"
label_CH4 = 'CH4'
tdms_file_path_CH4_new = 'ShortPulses/Alpha250/WorkAir#2_batch-0,5mbar_231118/Batch-0,5mbar_WorkAir#2(1)_CH4_231118.tdms'

"Foreground partial pressure fit parameter of main N2O line"
line_N2O = "Line 0Amplitude"
label_N2O = 'N2O'
tdms_file_path_N2O_new = 'ShortPulses/Alpha250/WorkAir#2_batch-0,5mbar_231118/Batch-0,5mbar_WorkAir#2(1)_N2O_231118.tdms'

"Foreground partial pressure fit parameter of main C12(CO2) line"
line_626 = "Line 0Amplitude"
label_626 = 'CO2 (626)'
tdms_file_path_626_new = 'ShortPulses/Alpha250/WorkAir#2_batch-0,5mbar_231118/Batch-0,5mbar_WorkAir#2(1)_CO2-4_231118.tdms'

"Foreground partial pressure fit parameter of main C13(CO2) line"
line_636 = "Line 0Amplitude"
label_636 = 'CO2 (636)'
tdms_file_path_636_new = 'ShortPulses/Alpha250/WorkAir#2_batch-0,5mbar_231118/Batch-0,5mbar_WorkAir#2(1)_CO2-3_231118.tdms'

label_ratio = 'd13C(CO2)'

"Background partial pressure fit parameter of selected line"
background = False #is there a background fittet at all?
line_bg1 = "Line 0Amplitude"
line_bg2 = "Line 1Amplitude"
label_bg1 = 'Background line 1'
label_bg2 = 'Background line 2'
tdms_file_path_bg_new = 'ShortPulses/Alpha250/WorkAir#2_batch-0,5mbar_231118/Batch-0,5mbar_WorkAir#2(1)_CO2-2_231118.tdms'

"concentration of main lines set in Hitran"
pCH4 = 1000 #[ppb]
pN2O = 100 #[ppb]
pCO2_fg = 100 #[ppm]
"Cell pressure set in Hitran [atm]"
cellPress = 0.005
"CO2 concentration of background lines set in Hitran [ppm]"
pCO2_bg = 4

"data clipping necesary?"
clip = False
sigma = 5

"""starting and ending data line to include"""
win_start = 0 #counted from beginning of file. First value = 0
win_end = 1 #counted from the back of file. Last value = 1

"""
Main Program
"""

"""load files if not yet loaded and make a safety copy."""
if 'tdms_file_path_CH4' not in globals() or tdms_file_path_CH4_new != tdms_file_path_CH4:
    tdms_file_path_CH4 = tdms_file_path_CH4_new
    tdms_file_CH4 = tdmsProcFuncs.readtdmsfile(tdms_file_path_CH4)
    tdms_file_CH4_copy = copy.deepcopy(tdms_file_CH4)

if 'tdms_file_path_N2O' not in globals() or tdms_file_path_N2O_new != tdms_file_path_N2O:
    tdms_file_path_N2O = tdms_file_path_N2O_new
    tdms_file_N2O = tdmsProcFuncs.readtdmsfile(tdms_file_path_N2O)
    tdms_file_N2O_copy = copy.deepcopy(tdms_file_N2O)

if 'tdms_file_path_626' not in globals() or tdms_file_path_626_new != tdms_file_path_626:
    tdms_file_path_626 = tdms_file_path_626_new
    tdms_file_626 = tdmsProcFuncs.readtdmsfile(tdms_file_path_626)
    tdms_file_CO2_626 = copy.deepcopy(tdms_file_626)

if 'tdms_file_path_636' not in globals() or tdms_file_path_636_new != tdms_file_path_636:
    tdms_file_path_636 = tdms_file_path_636_new
    tdms_file_636 = tdmsProcFuncs.readtdmsfile(tdms_file_path_636)
    tdms_file_CO2_636 = copy.deepcopy(tdms_file_636)

if 'tdms_file_path_bg' not in globals() or tdms_file_path_bg_new != tdms_file_path_bg:
    tdms_file_path_bg = tdms_file_path_bg_new
    tdms_file_bg = tdmsProcFuncs.readtdmsfile(tdms_file_path_bg)
    tdms_file_CO2_bg = copy.deepcopy(tdms_file_bg)

"""find group and channel names (i.e. "Fit Results" and "Line 1Amplitude, respectively)
via tdms_file.objects. Then access them as follows""" 

"""Time Stamp"""
timeStamp_fit_CH4_orig = tdms_file_CH4.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_CH4_orig = tdms_file_CH4.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

timeStamp_fit_N2O_orig = tdms_file_N2O.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_N2O_orig = tdms_file_N2O.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

timeStamp_fit_626_orig = tdms_file_626.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_626_orig = tdms_file_626.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

timeStamp_fit_636_orig = tdms_file_636.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_636_orig = tdms_file_636.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

timeStamp_fit_bg_orig = tdms_file_bg.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_bg_orig = tdms_file_bg.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

"""Find boarders to clip data files to the common time stamps"""
timeStamp_start = max(timeStamp_fit_CH4_orig[0],timeStamp_fit_N2O_orig[0],timeStamp_fit_626_orig[0],timeStamp_fit_636_orig[0],timeStamp_fit_bg_orig[0])
timeStamp_end = min(timeStamp_fit_CH4_orig[-1],timeStamp_fit_N2O_orig[-1],timeStamp_fit_626_orig[-1],timeStamp_fit_636_orig[-1],timeStamp_fit_bg_orig[-1])

ind_CH4 = [timeStamp_fit_CH4_orig.index(timeStamp_start)]
ind_CH4.append(timeStamp_fit_CH4_orig.index(timeStamp_end))
ind_N2O = [timeStamp_fit_N2O_orig.index(timeStamp_start)]
ind_N2O.append(timeStamp_fit_N2O_orig.index(timeStamp_end))
ind_626 = [timeStamp_fit_626_orig.index(timeStamp_start)]
ind_626.append(timeStamp_fit_626_orig.index(timeStamp_end))
ind_636 = [timeStamp_fit_636_orig.index(timeStamp_start)]
ind_636.append(timeStamp_fit_636_orig.index(timeStamp_end))
ind_bg = [timeStamp_fit_bg_orig.index(timeStamp_start)]
ind_bg.append(timeStamp_fit_bg_orig.index(timeStamp_end))

"""Transform TimeStamp into seconds after start"""
timeStamp_fit_CH4 = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_CH4_orig[ind_CH4[0]:ind_CH4[1]])
timeStamp_mon_CH4 = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_CH4_orig[ind_CH4[0]:ind_CH4[1]])
timeStamp_fit_N2O = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_N2O_orig[ind_N2O[0]:ind_N2O[1]])
timeStamp_mon_N2O = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_N2O_orig[ind_N2O[0]:ind_N2O[1]])
timeStamp_fit_626 = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_626_orig[ind_626[0]:ind_626[1]])
timeStamp_mon_626 = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_626_orig[ind_626[0]:ind_626[1]])
timeStamp_fit_636 = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_636_orig[ind_636[0]:ind_636[1]])
timeStamp_mon_636 = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_636_orig[ind_636[0]:ind_636[1]])
timeStamp_fit_bg = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_bg_orig[ind_bg[0]:ind_bg[1]])
timeStamp_mon_bg = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_bg_orig[ind_bg[0]:ind_bg[1]])

"""Check if time stamp of Monitoring Vals and Fit Results fit"""
if len(timeStamp_fit_626) != len(timeStamp_mon_626):
    warnings.warn('Time stamp of monitoring and fit values not same length!')
else:
    if not np.allclose(timeStamp_fit_626,timeStamp_mon_626, atol = 1, rtol = 0):
        warnings.warn('Time stamp of monitoring and fit values is shifted > 1sec.') 

"""Amplitudes of main lines"""
lineAmp_CH4 = tdms_file_CH4.object("Fit Results" , line_CH4).data[win_start:-win_end][ind_CH4[0]:ind_CH4[1]]
lineAmp_N2O = tdms_file_N2O.object("Fit Results" , line_N2O).data[win_start:-win_end][ind_N2O[0]:ind_N2O[1]]
lineAmp_626 = tdms_file_626.object("Fit Results" , line_626).data[win_start:-win_end][ind_626[0]:ind_626[1]]
lineAmp_636 = tdms_file_636.object("Fit Results" , line_636).data[win_start:-win_end][ind_636[0]:ind_636[1]]

"""Amplitudes of background lines"""
if background is True:
    lineAmp_bg1 = tdms_file_bg.object("Fit Results" , line_bg1).data[win_start:-win_end][ind_bg[0]:ind_bg[1]]
    lineAmp_bg2 = tdms_file_bg.object("Fit Results" , line_bg2).data[win_start:-win_end][ind_bg[0]:ind_bg[1]]
    """Calculate concentrations of background lines"""
    lineCons_bg1 = lineAmp_bg1*pCO2_bg
    lineCons_bg2 = lineAmp_bg2*pCO2_bg

"""Cell Pressure"""
p1Cell = tdms_file_626.object("Monitoring Vals" , "p1-Cell").data[win_start:-win_end][ind_626[0]:ind_626[1]]
#p1Cell = np.mean(p1Cell) #remove "noisy" pressure reading
#p1Cell = cellPress#for Zero Air evaluation

"""Other Monitoring Data"""
T1Cell = tdms_file_626.object("Monitoring Vals" , "T1-Cell").data[win_start:-win_end][ind_626[0]:ind_626[1]]
T_Cell_DAQVal = tdms_file_626.object("Monitoring Vals" , "T Cell DAQVal").data[win_start:-win_end][ind_626[0]:ind_626[1]]

"""Calculate main line ratio in per mil"""
line3to2_ratio = (1-(lineAmp_636/lineAmp_626)/(lineAmp_636[200]/lineAmp_626[200]))*1000

"""Calculate concentrations of main lines"""
lineCons_CH4 = lineAmp_CH4*pCH4*cellPress/p1Cell
lineCons_N2O = lineAmp_N2O*pN2O*cellPress/p1Cell
lineCons_626 = lineAmp_626*pCO2_fg*cellPress/p1Cell
lineCons_636 = lineAmp_636*pCO2_fg*cellPress/p1Cell

""" clip outliers from fitting errors if requested"""
if clip is True:
    flagratio = np.empty(5)
    flagratio.fill(np.nan)
    lineCons_626, flagratio[0] = tdmsProcFuncs.outlierFilt(lineCons_626,sigma)
    lineCons_636, flagratio[1] = tdmsProcFuncs.outlierFilt(lineCons_636,sigma)
    line3to2_ratio, flagratio[2] = tdmsProcFuncs.outlierFilt(line3to2_ratio,sigma)
    lineCons_CH4, flagratio[3] = tdmsProcFuncs.outlierFilt(lineCons_CH4,sigma)
    lineCons_N2O, flagratio[4] = tdmsProcFuncs.outlierFilt(lineCons_N2O,sigma)
    if background is True:
        lineCons_bg1, flagratio[3] = tdmsProcFuncs.outlierFilt(lineCons_bg1,sigma)
        lineCons_bg2, flagratio[4] = tdmsProcFuncs.outlierFilt(lineCons_bg2,sigma)
#    flagratio = np.nanmean(flagratio)

"""plot foreground data"""
#fig_fg, ax_cons, ax_rat = tdmsProcFuncs.plot3var_2yaxis(timeStamp_fit_626,lineCons_626,\
#                    lineCons_636,line3to2_ratio,label_626,label_636,label_ratio)
#ax_cons.set_ylabel('CO2 [ppm]')
#ax_cons.set_xlabel('time (sec)')
#ax_rat.set_ylabel('d13C(CO2) [per mil]')
#ttext = plt.title('CO2 parameters in cell')
#plt.setp(ttext, size='large', color='r', weight='bold')
#plt.text(3.0, 0.6, 'f(t) = exp(-t) sin(2 pi t)')

fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit_626,lineCons_626,'r',label=label_626)
plt.xlabel('Time [sec]')
plt.ylabel('CO2 [ppm]')
ax_1.legend(loc='lower left')
ax_1.grid(True)
ax_2 = ax_1.twinx()
ax_2.plot(timeStamp_fit_626,line3to2_ratio,'k',label=label_ratio)
plt.ylabel('d13C(CO2) [per mil]')
ax_2.legend(loc='lower right')
ttext = plt.title('CO2 concentration and d13C(CO2)')
plt.setp(ttext, size='large', color='r', weight='bold')

fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit_CH4,lineCons_CH4,'b',label=label_CH4)
plt.xlabel('Time [sec]')
plt.ylabel('CH4 [ppb]')
ax_1.legend(loc='lower left')
ax_1.grid(True)
ax_2 = ax_1.twinx()
ax_2.plot(timeStamp_fit_CH4,lineCons_N2O,'g',label=label_N2O)
plt.ylabel('N2O [ppb]')
ax_2.legend(loc='lower right')
ttext = plt.title('CH4 and N2O concentrations')
plt.setp(ttext, size='large', color='r', weight='bold')

"""plot background CO2 data"""
if background is True:
    fig_bg, ax_bg = tdmsProcFuncs.plot2var_1yaxis(timeStamp_fit_626,lineCons_bg1,\
             lineCons_bg2,label_bg1,label_bg2)
    ax_bg.set_ylabel('CO2 [ppm]')
    ax_bg.set_xlabel('time (sec)')
    ttext = plt.title('CO2 parameters of background')
    plt.setp(ttext, size='large', color='r', weight='bold')

"""Plot monitoring data"""
fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit_626,p1Cell,'r',label='Cell Pressure')
plt.xlabel('Time [sec]')
plt.ylabel('Cell Pressure [atm]')
ax_1.legend(loc='lower left')
ax_1.grid(True)
ax_2 = ax_1.twinx()
ax_2.plot(timeStamp_fit_626,T1Cell,'k',label='Cell Temperature')
plt.ylabel('Cell Temperature [°C]')
ax_2.legend(loc='lower right')
ttext = plt.title('Monitoring Data')
plt.setp(ttext, size='large', color='r', weight='bold')

#fig, ax_1 = plt.subplots()
#ax_1.plot(timeStamp_fit_626,T_Cell_DAQVal,'r',label='Raw Value DAQ')
#plt.xlabel('Time [sec]')
#plt.ylabel('Raw Value [V]')
#ax_1.legend(loc='lower left')
#ax_1.grid(True)
#ax_2 = ax_1.twinx()
#ax_2.plot(timeStamp_fit_626,T1Cell,'k',label='Cell Temperature')
#plt.ylabel('Cell Temperature [°C]')
#ax_2.legend(loc='lower right')
#ttext = plt.title('Cell temperature sensor: raw and transformed values')
#plt.setp(ttext, size='large', color='r', weight='bold')
#
#fig, ax_1 = plt.subplots()
#ax_1.plot(T1Cell,T_Cell_DAQVal,'r',label='Raw vs.transformed data')
#plt.xlabel('Cell Temperature [°C]')
#plt.ylabel('Raw Value [V]')
#ax_1.legend(loc='lower left')
#ax_1.grid(True)
#ttext = plt.title('Cell temperature sensor: raw vs. transformed values')
#plt.setp(ttext, size='large', color='r', weight='bold')
