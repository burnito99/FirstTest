# -*- coding: utf-8 -*-
"""
Script to to read CO2 fitting data and calculate the basic parameters based on 
the tdms file produced by the IRSENS Labview Software.
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import copy
import tdmsProcFuncs
from astropy.stats import sigma_clip

"""
input Data
"""
"Define file path with fiting results" 
tdms_file_path_new = 'ShortPulses/Redpitaya/StandAirBatch_070617/StandAirBatch_CO2_070617.tdms'

"Foreground partial pressure fit parameter of main C12(CO2) and C13(CO2) lines"
line_626 = "Line 0Amplitude"
line_636 = "Line 1Amplitude"
label_626 = 'CO2 (626)'
label_636 = 'CO2 (636)'
label_ratio = 'd13C(CO2)'

"Background partial pressure fit parameter of selected line"
background = False #is there a background fittet at all?
line_bg1 = "Line 10Amplitude"
line_bg2 = "Line 7Amplitude"
label_bg1 = 'Background line 1'
label_bg2 = 'Background line 2'

"CO2 concentration of main lines set in Hitran [ppm]"
pCO2_fg = 100
"Cell pressure set in Hitran [atm]"
cellPress = 0.005
"CO2 concentration of background lines set in Hitran [ppm]"
pCO2_bg = 4

"data clipping necesary?"
clip = False
sigma = 3

"""starting and ending data line to include"""
win_start = 0 #counted from beginning of file. First value = 0
win_end = 1 #counted from the back of file. Last value = 1

"other important parameters saved within 'Fit Results' (e.g. baseline parameters) "
parlist_co2 = [["P-Coeff 0"]]
parlist_co2.append(["P-Coeff 1"])
parlist_co2.append(["P-Coeff 2"])
parlist_co2.append(["P-Coeff 3"])
parlist_co2.append(["Zero Level"])

"""
Main Program
"""

"""load file if not yet loaded and make a safety copy."""
if 'tdms_file_path' not in globals() or tdms_file_path_new != tdms_file_path:
    tdms_file_path = tdms_file_path_new
    tdms_file = tdmsProcFuncs.readtdmsfile(tdms_file_path)
    tdms_file_CO2 = copy.deepcopy(tdms_file)

"""find group and channel names (i.e. "Fit Results" and "Line 1Amplitude, respectively)
via tdms_file.objects. Then access them as follows""" 

"""Time Stamp"""
timeStamp_fit_orig = tdms_file.object("Fit Results" , "t-stamp").data[win_start:-win_end]
timeStamp_mon_orig = tdms_file.object("Monitoring Vals" , "t-stamp").data[win_start:-win_end]

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
lineAmp_626 = tdms_file.object("Fit Results" , line_626).data[win_start:-win_end]
lineAmp_636 = tdms_file.object("Fit Results" , line_636).data[win_start:-win_end]

"""Amplitudes of background lines"""
if background is True:
    lineAmp_bg1 = tdms_file.object("Fit Results" , line_bg1).data[win_start:-win_end]
    lineAmp_bg2 = tdms_file.object("Fit Results" , line_bg2).data[win_start:-win_end]
    """Calculate concentrations of background lines"""
    lineCons_bg1 = lineAmp_bg1*pCO2_bg
    lineCons_bg2 = lineAmp_bg2*pCO2_bg

"""Cell Pressure"""
p1Cell = tdms_file.object("Monitoring Vals" , "p1-Cell").data[win_start:-win_end]
p1Cell = np.mean(p1Cell) #remove "noisy" pressure reading
p1Cell = cellPress#for Zero Air evaluation

"""Calculate main line ratio in per mil"""
line3to2_ratio = (1-(lineAmp_636/lineAmp_626)/(lineAmp_636[200]/lineAmp_626[200]))*1000

"""Calculate concentrations of main lines"""
lineCons_626 = lineAmp_626*pCO2_fg*cellPress/p1Cell
lineCons_636 = lineAmp_636*pCO2_fg*cellPress/p1Cell

""" clip outliers from fitting errors if requested"""
if clip is True:
    flagratio = np.empty(5)
    flagratio.fill(np.nan)
    lineCons_626, flagratio[0] = tdmsProcFuncs.outlierFilt(lineCons_626,sigma)
    lineCons_636, flagratio[1] = tdmsProcFuncs.outlierFilt(lineCons_636,sigma)
    line3to2_ratio, flagratio[2] = tdmsProcFuncs.outlierFilt(line3to2_ratio,sigma)
    if background is True:
        lineCons_bg1, flagratio[3] = tdmsProcFuncs.outlierFilt(lineCons_bg1,sigma)
        lineCons_bg2, flagratio[4] = tdmsProcFuncs.outlierFilt(lineCons_bg2,sigma)
    flagratio = np.nanmean(flagratio)

"""Unpack other parameters"""
parvalues_co2 = []
for x in range(len(parlist_co2)):
    parvalues_co2.append(tdms_file.object("Fit Results" , *parlist_co2[x]).data[win_start:-win_end])

"""plot foreground CO2 data"""
fig_fg, ax_cons, ax_rat = tdmsProcFuncs.plot3var_2yaxis(timeStamp_fit,lineCons_626,\
                    lineCons_636,line3to2_ratio,label_626,label_636,label_ratio)
ax_cons.set_ylabel('CO2 [ppm]')
ax_cons.set_xlabel('time (sec)')
ax_rat.set_ylabel('d13C(CO2) [per mil]')
ttext = plt.title('CO2 parameters in cell')
plt.setp(ttext, size='large', color='r', weight='bold')
#plt.text(3.0, 0.6, 'f(t) = exp(-t) sin(2 pi t)')

fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit,lineCons_626,'r',label=label_626)
plt.xlabel('Time [sec]')
plt.ylabel('CO2 [ppm]')
ax_1.legend(loc='lower left')
ax_1.grid(True)
ax_2 = ax_1.twinx()
ax_2.plot(timeStamp_fit,line3to2_ratio,'k',label=label_ratio)
plt.ylabel('d13C(CO2) [per mil]')
ax_2.legend(loc='lower right')
ttext = plt.title('CO2 Concentrations in Batch Mode')
plt.setp(ttext, size='large', color='r', weight='bold')

"""plot background CO2 data"""
if background is True:
    fig_bg, ax_bg = tdmsProcFuncs.plot2var_1yaxis(timeStamp_fit,lineCons_bg1,\
             lineCons_bg2,label_bg1,label_bg2)
    ax_bg.set_ylabel('CO2 [ppm]')
    ax_bg.set_xlabel('time (sec)')
    ttext = plt.title('CO2 parameters of background')
    plt.setp(ttext, size='large', color='r', weight='bold')
        
"""plot other data"""
fig_pol, ax_0, ax_1 = tdmsProcFuncs.plot3var_2yaxis(timeStamp_fit,parvalues_co2[0],\
                    parvalues_co2[1],parvalues_co2[2],parlist_co2[0],\
                    parlist_co2[1],parlist_co2[2])
ax_0.set_xlabel('time (sec)')
ax_0.set_ylabel('y1 and y2 values')
ax_1.set_ylabel('y3 and y4 values')
ttext = plt.title('Baseline Polynom Parameters of 4.3 um laser')
plt.setp(ttext, size='large', color='r', weight='bold')
