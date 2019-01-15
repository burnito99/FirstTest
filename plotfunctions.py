# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:47:31 2017

@author: scsi
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.dates as mdates
import allantools as at
from scipy.fftpack import fft, fftfreq, fftshift


#def xy_plot(x,y,label):
#    plt.figure()
#    fig, ax = plt.subplots()
#    ax.plot(x,y,label=label)
#    ax.legend(loc='best')
#    ax.grid(True)
#    plt.show()
"""
x = 1D array y=nD array, 
label: label of the the lines, 
plotaxis (n,m) array with the lines corresponding

"""
def shape_array(A):
    if np.size(A) == 0:
        size = 0
    else:
        size = np.shape(A)[0]
    return size


def ma_concatination(A,B):
    size_A = shape_array(A)
    size_B = shape_array(B)
    y = [None]*(int(size_A) + int(size_B))
    for i in range(size_A):
        y[i] = A[i]
    for i in range(size_B):
        y[size_A + i] = B[i]
    return y

def check_mask(x):
    if True in np.ma.getmaskarray(x):
        xx = x[~x.mask]
    else:
        xx = x
    return xx

def ratio_builder_allan(y,pl, labels):
    data = []
    label = []
    k = 1
    if pl == []:
        pass
    else:
        
        for task in pl:
            if ":" in task:
                A = task.split(":")
                ratio = 1000.0*(1.0 - np.array(y[int(A[0])])*np.mean(np.array(y[int(A[1])]))/(np.array(y[int(A[1])])* np.mean(np.array(y[int(A[0])]))))
#                ratio = np.ma.array(ratio)
#                ratio = ratio[~ratio.mask]
                data.append(ratio)
                name = 'Ratio ' + str(k)
                label.append(name)
                k += 1
            else:
                data.append(y[int(task)])
                label.append(labels[int(task)])
    return data, label      
        
def ratio_builder(y,plotaxis, legend,k):
    B = []
    Ratio = []
    for task in plotaxis:
        yy = int(np.shape(y)[0])
        if ":" in task:
            A = task.split(":")
            ratio = 1000.0*(1.0 - y[int(A[0])]*np.mean(y[int(A[1])])/(y[int(A[1])]* np.mean(y[int(A[0])])))
            Ratio.append(ratio)
            B.append(str(yy))
            name = "Ratio " + str(k)
            legend.insert(yy,name)
#            legend.append(name)
            yy += 1
            k += 1
        else:
            B.append(task)
    return Ratio, B, legend, k

def ratio_generator(y,y2,plotaxis, legend):
    k = 1
    
    ratio, B, legend, k = ratio_builder(y, plotaxis[0][0],legend, k)
    if ratio == []:
        pass
    else:
        if np.size(y[0]) == np.size(ratio[0]):
            y = ma_concatination(y,ratio)
        else:
            pass
    
    plotaxis[0][0] = B
    
    ratio, B, legend, k = ratio_builder(y2, plotaxis[0][1], legend, k)
    if ratio == []:
        pass
    else:
        if np.size(y2[0]) == np.size(ratio[0]):
            y = ma_concatination(y,ratio)
        else:
            pass
    
    plotaxis[0][1] = B
    
    ratio, B, legend, k = ratio_builder(y, plotaxis[1][0], legend, k)
    if ratio == []:
        pass
    else:
        if np.size(y[0]) == np.size(ratio[0]):
            y = ma_concatination(y,ratio)
        else:
            pass
    
    plotaxis[1][0] = B
    
    ratio, B, legend, k = ratio_builder(y2, plotaxis[1][1], legend, k)
    if ratio == []:
        pass
    else:
        if np.size(y2[0]) == np.size(ratio[0]):
            y = ma_concatination(y,ratio)
        else:
            pass
    
    plotaxis[1][1] = B
    return y, y2, plotaxis, legend

def colormap_generator(plotaxis):
    plot_number = np.size(plotaxis) + 4
    color = iter(plt.cm.rainbow(np.linspace(0,1,plot_number)))
    return color
    
def VoltsToTempConvert(voltsdata):
    """converts volts measured at a temperature sensor into temperature"""
    temp = 1/(1.1279*10**(-3)+2.3429*10**(-4)*np.log(voltsdata/0.0001)\
              +8.7298*10**(-8)*(np.log(voltsdata/0.0001))**3)-273.15
    return temp                     
            
        
    
    
def xy_plot(x,y,x2,y2, l_s, m, label, title, x_lab, yl_lab,yr_lab, plotaxis = np.array([[["0"],[""]],[[""],["1"]]]), twoyaxis = False, legends = True):
    y, y2, plotaxis, label = ratio_generator(y,y2,plotaxis, label)
    colors = colormap_generator(plotaxis)
    if np.shape(y)[1] != 0:
        y_size = np.shape(y)[0]
    else:
        y_size = 0
    if np.shape(y2)[1] != 0:
        y2_size = np.shape(y2)[0]
    else:
        y2_size = 0
    if twoyaxis == False:      
        fig, ax = plt.subplots()
        ax.set_xlabel(x_lab)
        ax.set_ylabel(yl_lab)
        plt.title(title)
        #plot the given axis
        if "" in plotaxis[0][0]:
            pass
        else:
            for i in range(np.size(plotaxis[0][0])):
                col = next(colors)
                ax.plot(x,y[int(plotaxis[0][0][i])], c = col, ls = l_s[0][i%np.size(l_s[0])], marker = m[0][i%np.size(m[0])], label=label[int(plotaxis[0][0][i])]
                )
        if "" in plotaxis[0][1]:
            pass
        else:
            for i in range(np.size(plotaxis[0][1])):
                col = next(colors)
                ax.plot(x2,y2[int(plotaxis[0][1][i])], c = col, ls = l_s[0][i%np.size(l_s[0])], marker = m[0][i%np.size(m[0])], label=label[y_size + int(plotaxis[0][1][i])]
                )
        if legends == True:            
            ax.legend(loc='best')
        else:
            pass
        ax.grid(True)


    else:
        fig, ax_l = plt.subplots()
        if "" in plotaxis[0][0]:
            pass
        else:
            for i in range(np.size(plotaxis[0][0])):
                col = next(colors)
                ax_l.plot(x,y[int(plotaxis[0][0][i])], c = col, ls = l_s[0][i%np.size(l_s[0])], marker = m[0][i%np.size(m[0])], label=label[int(plotaxis[0][0][i])]
                )
        if "" in plotaxis[0][1]:
            pass
        else:
            for i in range(np.size(plotaxis[0][1])):
                col = next(colors)
                ax_l.plot(x2,y2[int(plotaxis[0][1][i])], c = col, ls = l_s[0][i%np.size(l_s[0])], marker = m[0][i%np.size(m[0])], label=label[y_size + int(plotaxis[0][1][i])]
                )
            
        ax_r = ax_l.twinx()
        if "" in plotaxis[1][0]:
            pass
        else:
            for i in range(np.size(plotaxis[1][0])):
                col = next(colors)
                ax_r.plot(x,y[int(plotaxis[1][0][i])], c = col, ls = l_s[1][i%np.size(l_s[1])], marker = m[1][i%np.size(m[1])], label=label[int(plotaxis[1][0][i])]
                )
        if "" in plotaxis[1][1]:
            pass
        else:
            for i in range(np.size(plotaxis[1][1])):
                col = next(colors)
                ax_r.plot(x2,y2[int(plotaxis[1][1][i])], c = col, ls = l_s[1][i%np.size(l_s[1])], marker = m[1][i%np.size(m[1])], label=label[y_size +int(plotaxis[1][1][i])]
                )
        ax_l.set_xlabel(x_lab)
        ax_l.set_ylabel(yl_lab)
        ax_r.set_ylabel(yr_lab)
        ax_l.grid(True)
#        ax_r.grid(True)
        if legends == True:
            lines_l, labels_l = ax_l.get_legend_handles_labels()
            lines_r, labels_r = ax_r.get_legend_handles_labels()
            ax_r.legend(lines_l + lines_r, labels_l + labels_r, loc = 'best'#,bbox_to_anchor = (1,1)
                        )
        plt.title(title)
        plt.show()



#def oallan_plot(y, plotaxis, rate, taus, style, marker, x_lab, yl_lab, label, title, legends = True):
#    y = ratio_builder_allan(y, plotaxis)
##    colors = colormap_generator(plotaxis)
#    if np.shape(y)[1] != 0:
#        y_size = np.shape(y)[0]
#    else:
#        y_size = 0
#
#    plt.figure()
#    plt.subplot(111, xscale="log", yscale="log", xlabel = x_lab, ylabel = yl_lab, title = title)
#    for i in range(y_size): 
##        col = next(colors)
#        (t2, ad, ade, adn) = at.oadev(y[i], rate=rate, data_type="freq", taus='all')
#        plt.loglog(t2, ad,"g.",# c = col,# ls = str(style), marker = "o",
#                   label = label[i])
#    if legends == True:
#        plt.legend(loc = 'best')
#    plt.show()
def oallan_plot(y, plotaxis, rate, taus, style, marker, x_lab, yl_lab, label, title, legends = True):
    y, label = ratio_builder_allan(y, plotaxis, label)
    colors = colormap_generator(plotaxis)
    if np.size(y) != 0:
        y_size = np.shape(y)[0]
    else:
        y_size = 0

    plt.figure()
    ax = plt.subplot(111, xscale="log", yscale="log", xlabel = x_lab, ylabel = yl_lab, title = title)
    for i in range(y_size):
        yy = check_mask(y[i])
        col = next(colors)
        (t2, ad, ade, adn) = at.oadev(yy, rate=rate, data_type="freq", taus='all')
        mini = np.argmin(ad[:int(np.size(ad)/2)])
        allan_minimum = '{:1.1e}'.format(ad[mini])
        time_minimum = '{:1.1e}'.format(t2[mini])
        label[i] = label[i].replace("Amplitude","") + " ["+ allan_minimum +"; "+ time_minimum + "]"
        
        ax.loglog(t2, ad,"g.",c = col, ls = style[i%np.size(style)], marker = marker[i%np.size(marker)],
                   label = label[i])#label[i])
    if legends == True:
        ax.legend(loc = 'best')
    plt.show()

    
def moving_allan(y,x, length, size, rate, plotaxis,yl_lab,title, label):
    y_new, label = ratio_builder_allan(y, plotaxis, label)

    if np.shape(y_new)[1] != 0:
        y_size = np.shape(y_new)[0]
    else:
        y_size = 0

    k = int((np.shape(y_new)[1]- length)/size) + 1
    x_allan =(length/2.0 + np.arange(0,k)*size)/rate + x[0]
    minimas = np.zeros((2,y_size,k))
    for i in range(y_size):
        yy = check_mask(y_new[i])
        for j in range(k):
            (t2, ad, ade, adn) = at.oadev(yy[j*size:length + j*size], rate=rate, data_type="freq", taus='all')
            minimas[0][i][j] = (t2[np.argmin(ad[:int(length/2)])])
            minimas[1][i][j] = (np.min(ad[:int(length/2)]))
        H, xedges, yedges = (np.histogram2d(minimas[1][i],minimas[0][i],bins = 10))
        H /= k
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        X, Y = np.meshgrid(xedges[:-1]+ 0.1*xedges[0],yedges[:-1]+0.1*yedges[0])
        X = X.flatten('F')
        Y = Y.flatten('F')
        H = H.flatten()
        Z = np.zeros_like(X)
        dx = 0.8*abs(xedges[1]- xedges[0])
        dy = 0.8*abs(yedges[1]- yedges[0])
        ax.bar3d(X,Y,Z,dx,dy,H,zsort = 'average')
        ax.set_xlabel('Allan minimum')
        ax.set_ylabel('Integration time (sec)')
        ax.set_zlabel('normalized frequency')
        fig.show()
        plt.figure()   
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212, sharex = ax1)
        ax3 = ax2.twinx()
        x = np.linspace(0,np.size(yy)/rate,num = np.size(yy))
        ax1.plot(x,yy)
#        ax1.set_xlim
        ax2.set_yscale('log')
        ax2.errorbar(x_allan, minimas[1][0],xerr = length/(2*rate), c = 'r',ls = '', label = 'Allan minima')
#        ax2.set_xlim([min(x),max(x)])
        ax3.plot(x_allan, minimas[0][0], c = 'b', ls = '', marker = 'o', label = "Integration time")
        ax1.set_ylabel(yl_lab)
        ax3.set_ylabel("Time (sec)")
        ax2.set_ylabel("Allan minimum")
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines3, labels3 = ax3.get_legend_handles_labels()
        ax3.legend(lines2 + lines3, labels2 + labels3, loc = 4)
        plt.title(title)
        plt.show()
    
    





def TimeStampTransform(timeStamp):
    """Transform timeStamp into seconds after start"""
    timeStamp = mdates.date2num(timeStamp)
    timeStamp = timeStamp-timeStamp[1] #This is time expressed in days
    timeStamp = timeStamp*3600*24 #Convert to Seconds 
    return timeStamp



def powerspectra(time,y, apodisation, zerofilling, yl_lab, yr_lab, title):#, plotaxis, label):
   # y, label = ratio_builder_allan(y, plotaxis, label)
    for i in range(np.shape(y)[0]):
        amplitude = y[i]
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

        
    
    
def autocor(time, y,yl_lab, x_lab, title, plotaxis,label):
    y, label = ratio_builder_allan(y, plotaxis, label)
    for i in range(np.shape(y)[0]):
        y1 = check_mask(y[i])    
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

    