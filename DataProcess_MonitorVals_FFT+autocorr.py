# -*- coding: utf-8 -*-
"""
Script to to read and plot monitoring parameters from a tdms file produced by 
the IRSENS Labview Software.
"""

import matplotlib.pyplot as plt
import numpy as np
import tdmsProcFuncs
import copy
import plotfunctions
from scipy.fftpack import fft, fftfreq

def check_mask(x):
    if True in np.ma.getmaskarray(x):
        xx = x[~x.mask]
    else:
        xx = x
    return xx

"""plot power spectra"""
def powerspectra(time,y, apodisation, zerofilling, yl_lab, yr_lab, title):
    amplitude = y
    length = amplitude.size
    duration = np.ptp(time)
    timestep = duration/length
    duration *= (1+zerofilling)
    apod = np.exp(-apodisation*np.linspace(0,duration, num = int(amplitude.size))/duration)
    amplitude = (amplitude-np.mean(amplitude))*apod
    amplitude = np.ma.concatenate((amplitude, np.zeros(int(amplitude.size*zerofilling))))
    
    length = int(amplitude.size)
    frequency = fftfreq(length ,d = timestep)
    fftampl = fft(amplitude)
    #        if length/2.0 is int:
    freq = frequency[1:int(length/2)]
    amp = fftampl[1:int(length/2)]
    #        else:
    #            freq = frequency[1:(length-1)/2]
    #            amp = fftampl[1:(length-1)/2]
    
    powerspec = np.absolute(amp)**2
    
    fig, (ax1, ax2) = plt.subplots(2,1)
    time = np.linspace(0,duration,num = length)
    plt.title(title)
    ax1.plot(time,amplitude)
    ax1.set_xlabel("Time (sec)")
    ax1.set_ylabel(yl_lab)
    ax2.semilogy(freq, powerspec)
    ax2.set_ylabel(yr_lab)
    ax2.set_xlabel("Frequency (Hz)")
    
    fig.show()
    
"plot auto correlation"
def autocor(time, y,yl_lab, x_lab, title):
    y1 = check_mask(y)    
    time -= np.min(time)
    y1 -= np.mean(y1)
    cor = np.correlate(y1,y1, mode = "full")
    #    x = np.linspace(- cor.size/2, cor.size/2, num = cor.size)
    autocorol = cor[int((cor.size-1)/2):]/cor[int((cor.size-1)/2)]
    x = np.linspace(0, 1, num = autocorol.size)*max(time)
        
    fig, ax = plt.subplots()
    ax.plot(x,autocorol)
    ax.set_ylabel(yl_lab)
    ax.set_xlabel(x_lab)
    plt.title(title)
    fig.show()


"""
input Data
"""

"Define file path with fiting results" 
tdms_file_path_new = 'ShortPulses/Alpha250/ADC_DataFromPhilip/ADC_fix.tdms'

"""Select paremeter names to plot. Second value indicates whether the data has
to be transformed from Volt to °C or not."""

parlist = [[["unnamed 15 DAQVal"],[False]]]#0
#parlist = [[["p1-Cell"],[False]]]#0
parlist.append([["unnamed 16 DAQVal"],[False]])#1
parlist.append([["unnamed 19 DAQVal"],[False]])#2
parlist.append([["unnamed 20 DAQVal"],[False]])#3
#parlist.append([["T Box DAQVal"],[True]])#4
#parlist.append([["Vici Valve Pos DAQVal"],[False]])#5
#parlist.append([["p Turbo DAQVal"],[False]])#6
#parlist.append([["Cell Entry Valve DAQVal"],[False]])#7
#parlist.append([["ByPass Valve DAQVal"],[False]])#8
#parlist.append([["Turbo Valve DAQVal"],[False]])#9

par_ind = 0
ind_start = 0
ind_end = -1

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
timeStamp_mon_orig = tdms_file.object("Monitoring Values" , "t-stamp").data
#timeStamp_fit_orig = tdms_file.object("Fit Results" , "t-stamp").data

"""Transform TimeStamp into seconds after start"""
timeStamp_mon = tdmsProcFuncs.TimeStampTransform(timeStamp_mon_orig)
#timeStamp_fit = tdmsProcFuncs.TimeStampTransform(timeStamp_fit_orig)

"""Unpack parameter values and transform to temperature where necessary"""
parvalues = []
for x in range(len(parlist)):
    parvalues.append(tdms_file.object("Monitoring Values" , *parlist[x][0]).data)
    if parlist[x][1] == [True]:
        parvalues[x] = tdmsProcFuncs.VoltsToTempConvert(parvalues[x])

y = parvalues[par_ind][ind_start:ind_end]
time = timeStamp_mon[ind_start:ind_end]
zerofilling = 0
apodisation = 0
title_ps = 'powerspectrum'
yl_lab_ps = 'signal anomaly'
yr_lab = 'amplitude'
title_ac = 'autocorrelation'
x_lab = 'time [sec]'
yl_lab_ac = 'correlation factor'

"overview plot"
plt.figure()
plt.plot(time,y)
plt.ylabel('signal')
plt.xlabel('Time [sec]')

powerspectra(time,y, apodisation, zerofilling, yl_lab_ps, yr_lab, title_ps)

#autocor(time, y,yl_lab_ac, x_lab, title_ac)

