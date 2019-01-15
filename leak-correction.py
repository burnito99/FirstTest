# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 15:45:50 2018

@author: beb
"""

"concentrations in gas that leaks into sample"
pCH4_leak = 25000 #[ppb]
pN2O_leak = 3500 #[ppb]
pCO2_leak = 3430 #[ppm]

"transparancy for uncorrected data"
Alpha = 0.5

lineCons_626_corr = lineCons_626-(p1Cell-p1Cell[0])/p1Cell[0]*pCO2_leak
lineCons_CH4_corr = lineCons_CH4-(p1Cell-p1Cell[0])/p1Cell[0]*pCH4_leak
lineCons_N2O_corr = lineCons_N2O-(p1Cell-p1Cell[0])/p1Cell[0]*pN2O_leak

"CO2 plot"
fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit_626,lineCons_626,'r',alpha=Alpha,label='CO2 raw')
ax_1.plot(timeStamp_fit_626,lineCons_626_corr,'r',label='CO2 corrected')
plt.xlabel('Time [sec]')
plt.ylabel('CO2 [ppm]')
ax_1.legend(loc='best')
ax_1.grid(True)
ttext = plt.title('CO2 raw and leak-corrected concentration')
plt.setp(ttext, size='large', color='r', weight='bold')

"CH4 and N2O plot"
fig, ax_1 = plt.subplots()
ax_1.plot(timeStamp_fit_CH4,lineCons_CH4,'b',alpha=Alpha,label='CH4 raw')
ax_1.plot(timeStamp_fit_CH4,lineCons_CH4_corr,'b',label='CH4 corrected')
plt.xlabel('Time [sec]')
plt.ylabel('CH4 [ppb]')
ax_1.legend(loc='lower right')
ax_1.grid(True)
ax_2 = ax_1.twinx()
ax_2.plot(timeStamp_fit_CH4,lineCons_N2O,'g',alpha=Alpha,label='N2O raw')
ax_2.plot(timeStamp_fit_CH4,lineCons_N2O_corr,'g',label='N2O corrected')
plt.ylabel('N2O [ppb]')
ax_2.legend(loc='upper left')
ttext = plt.title('CH4 and N2O concentrations (raw and leak-corrected)')
plt.setp(ttext, size='large', color='r', weight='bold')
