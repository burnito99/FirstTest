# -*- coding: utf-8 -*-

"""
Script to calculate Allan Deviations and create Allan Plots
 Based on noise_color_demo.py by Anders Wallin (anders.e.e.wallin "at" gmail.com)
 v1.0 2014 January
"""

import matplotlib.pyplot as plt
import numpy as np
import allantools
import tdmsProcFuncs

"""
Inputs
"""

"Set Time Stamp"
ind_min = 0
ind_max = -1

#timeStamp_win = timeStamp_fit[ind_min:ind_max]
#timeStamp_win = timeStamp_fit[ind_min:ind_max]
timeStamp_win = timeStamp_fit_626[ind_min:ind_max]

"Set Time Series"
timeseries_win = lineCons_CH4[ind_min:ind_max]
#timeseries_win = lineCons_626[ind_min:ind_max]
#timeseries_win = line3to2_ratio[ind_min:ind_max]

"maximum of Allan devition curve"
Allan_max = 200 #[sec]

"Moving Allan Deviation Analysis"
allanMove = False #Do you want to derive Allan-minimum with a moving data window?
allanWinSize = 3000 #Window Size in seconds
winStep = 1/2 #step size from one window to the next

"Plot reference noise slopes"
plotnoise = False

"""
Main
"""

"Derive average sampling frequency"
delta_time = np.diff(timeStamp_win)
sampRate_av = 1/delta_time.mean()

"create figure with errorbars"
plt.figure()
ax = plt.axes()
ax.set_xscale("log")
ax.set_yscale("log")

#allanPar_adev = tdmsProcFuncs.plotallan_adev(plt, timeseries_win, sampRate_av, 'all', 'y.') #classic Allan deviation
#allanPar_oadev = tdmsProcFuncs.plotallan_oadev(plt, timeseries_win, sampRate_av, 'all', 'g') #considered as state-of-the-art
#allanPar_mdev = tdmsProcFuncs.plotallan_mdev(plt, timeseries_win, sampRate_av, 'all', 'm.') #able to distinguish between white and flicker noise

(t2, ad, ade, adn) = allantools.oadev(timeseries_win, rate=sampRate_av, data_type="freq", taus='all')#considered as state-of-the-art Allan Deviation algorithm
t2_ind_max = tdmsProcFuncs.bisection(t2,Allan_max)
ad_cut = ad[:t2_ind_max]
t2_cut = t2[:t2_ind_max]
ade_cut = ade[:t2_ind_max]
ax.errorbar(t2_cut, ad_cut, yerr = ade_cut, fmt = 'g')
plt.xlim(1, Allan_max)
ad_min_y = round(ad_cut.min(),4)#save Allan minimum
ad_min_x = round(t2_cut[np.argmin(ad_cut)])#save corresponding averaging time/tau
string = "Allan Minimum " + str(ad_min_y) + "\nIntegration time " + str(ad_min_x)
plt.text(ad_min_x, ad_min_y+(ad_cut.max()-ad_min_y)/2, string)
plt.ylabel('Allan Deviation [per mil]')
plt.xlabel('Time [sec]')
ttext = plt.title('Allan Deviation')
plt.setp(ttext, size='large', color='r', weight='bold')

"plot noise with known characteristics"
if plotnoise is True:
    noise_white = allantools.noise.white(len(timeStamp_win),allanPar_oadev[1][0]**2) #generate white noise data with slightly lower noise density as analyzed data
    allanPar_oadev_noise = tdmsProcFuncs.plotallan_oadev(plt, noise_white, sampRate_av, 'all', 'k.') #plot white noise data for reference
    tdmsProcFuncs.plotline(plt, -0.5, allanPar_oadev[1][0], allanPar_oadev[0], 'r')#plot -1/2 line representing white noise characteristic
# Dito as above but for pink and random walk noise
#    noise_pink = allantools.noise.pink(len(timeStamp_win))
#    noise_pink[:] = [x*allanPar_oadev[1][0] for x in noise_pink]
#    allanPar_oadev_noise = tdmsProcFuncs.plotallan_oadev(plt, noise_pink, sampRate_av, 'all', 'm.') #plot pink noise for reference
#    noise_randWalk = allantools.noise.brown(len(timeStamp_win),allanPar_oadev[1][0]**2)
#    allanPar_oadev_noise = tdmsProcFuncs.plotallan_oadev(plt, noise_randWalk, sampRate_av, 'all', 'y.') #plot random walk noise for reference

if allanMove is True:
    "get number of windows to fit into data set"
    winNum = int(1/winStep*(timeStamp_win[-1]-timeStamp_win[0]-allanWinSize/2)/allanWinSize)#number of windows that fit in data set
    if winNum == 0:
        warnings.warn('Allan Window does not fit into given data set')
    else:
        res = np.zeros([winNum, 4])#array to store window middle point, Allan minimum, corresponding tau, Allan start value
        timeStamp_win_win = timeStamp_win-timeStamp_win[0]
        "treat special case first window"
        ind0 = 0 #starting index
        ind1 = tdmsProcFuncs.bisection(timeStamp_win_win,allanWinSize/2)#find middle index
        ind2 = tdmsProcFuncs.bisection(timeStamp_win_win,allanWinSize)#find end index
        res[0,0] = timeStamp_win[ind1]#save middle time stamp
        "calculate allan minimum and its tau"
        delta_time = np.diff(timeStamp_win[ind0:ind2])
        sampRate_av = 1/delta_time.mean()
        (t, oad, oade, oadn) = allantools.oadev(timeseries_win[ind0:ind2], rate=sampRate_av, data_type="freq", taus='all')
        oad = oad[0:int(len(oad)/2)] #ingnore "uncertain" second half of allan deviation
        t = t[0:len(oad)]
        res[0,1] = oad.min()#save Allan minimum
        res[0,2] = t[np.argmin(oad)]#save corresponding averaging time/tau
        res[0,3] = oad[0] #save Allan Start value
        
        "loop over other windows"
        for xx in range(1,winNum): 
            winStart = xx*winStep*allanWinSize
            ind0 = tdmsProcFuncs.bisection(timeStamp_win_win,winStart)#find start index
            ind1 = tdmsProcFuncs.bisection(timeStamp_win_win,winStart+allanWinSize/2)#find middle index
            ind2 = tdmsProcFuncs.bisection(timeStamp_win_win,winStart+allanWinSize)#find end index
            res[xx,0] = timeStamp_win[ind1]#save middle time stamp
            delta_time = np.diff(timeStamp_win[ind0:ind2])
            sampRate_av = 1/delta_time.mean()
            (t, oad, oade, oadn) = allantools.oadev(timeseries_win[ind0:ind2], rate=sampRate_av, data_type="freq", taus='all')
            oad = oad[0:int(len(oad)/2)] #ingnore "uncertain" second half of allan deviation 
            t = t[0:len(oad)]
            res[xx,1] = oad.min()#save Allan minimum
            res[xx,2] = t[np.argmin(oad)]#save corresponding averaging time/tau
            res[xx,3] = oad[0] #save Allan Start value
            
        "make plot of Allan minimum and its tau based on the moving allan windows"
        fig_movA, ax_minA, ax_tA = tdmsProcFuncs.plot3var_2yaxis(res[:,0],\
                    res[:,1],res[:,3],res[:,2],'Allan-Minimum','Allan-Start','Integration Time')
        ax_minA.set_ylabel('Allan Deviation [per mil]')
        ax_minA.set_yscale('log')
        ax_minA.set_xlabel('time (sec)')
        ax_tA.set_ylabel('integration time [sec]')
        ttext = plt.title('Allan-Minimum and the Corresponding Integration Time')
        plt.setp(ttext, size='large', color='r', weight='bold')

""" Example how it works in "phase" world:
plt.figure()
plt.subplot(111, xscale="log", yscale="log")

timeseries_sum = np.cumsum(timeseries_win) # integrate to get  phase

tdmsProcFuncs.plotallan_oadev_phase(plt, timeseries_sum, delta_time_av, 'all', 'go')
tdmsProcFuncs.plotallan_adev_phase(plt, timeseries_sum, delta_time_av, 'all', 'yo')
tdmsProcFuncs.plotallan_mdev_phase(plt, timeseries_sum, delta_time_av, 'all', 'mo')
"""


""" Stuff from Demo File:
plt.figure()
t = [tt for tt in np.logspace(0, 3, 50)]  # tau values from 1 to 1000
plt.subplot(111, xscale="log", yscale="log")
N = 100000
# Colors: http://en.wikipedia.org/wiki/Colors_of_noise

# Brownian a.k.a random walk  frequency => sqrt(tau) ADEV
print("Random Walk frequency noise - should have sqrt(tau) ADEV")
freq_rw = noise.brown(N)
phase_rw_rw = np.cumsum(noise.brown(N))  # integrate to get  phase
plotallan_mdev(plt, freq_rw, 1, t, 'm.')
plotallan_mdev_phase(plt, phase_rw_rw, 1, t, 'mo',label='random walk frequency')
plotline(plt, +0.5, t, 'm',label="f^(+1/2)")

plotallan_adev(plt, freq_rw, 1, t, 'y.')
plotallan_adev_phase(plt, phase_rw_rw, 1, t, 'yo')

plotallan_oadev(plt, freq_rw, 1, t, 'g.')
plotallan_oadev_phase(plt, phase_rw_rw, 1, t, 'go')
"""