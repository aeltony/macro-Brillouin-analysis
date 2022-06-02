#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:07:33 2021

@author: ame38
"""
import os
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
plt.style.use('/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Data analysis/Brillouin human analysis/myPlotStyle.mplstyle')
plt.rcParams['figure.max_open_warning'] = False
plt.rcParams['pdf.fonttype'] = 42

#### Plot raw spectra images from Andor camera
def plotAndorImage(path, rawData, exp, scan, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'Spectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'Spectra' + '/')
    frames = rawData[exp][scan]['attrs']['paramtree/Scan/Frame Number']
    for i in range(frames):
        image = rawData[exp][scan]['AndorImage'][i]
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + ': frame %d' %i)
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest')
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view

        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'Spectra' + '/' \
                + exp + '_' + scan + r'_frame_%d_raw.png' %i, \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot calibration spectra images from Andor camera
def plotCalAndorImage(path, rawData, exp, scan, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'CalSpectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'CalSpectra' + '/')
    frames = rawData[exp][scan]['attrs']['paramtree/Scan/Frame Number']
    calFrames = rawData[exp][scan]['CalFreq'].shape[0]
    for l in range(calFrames):
        image = rawData[exp][scan]['AndorImage'][frames + l]
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + ': calFrame %d \n EOM setting: %.3f GHz' \
            %(l, rawData[exp][scan]['CalFreq'][l]))
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest')
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view

        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'CalSpectra' + '/' \
                + exp + '_' + scan + r'_calFrame_%d_raw.png' %l, \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot Brillouin spectrograph (binned)
def plotSpectrograph(path, rawData, procData, exp, scan, saveFig=False):
    frames = rawData[exp][scan]['attrs']['paramtree/Scan/Frame Number']
    pixels = procData[exp][scan]['RawSpecList'].shape[1]
    scanRange = rawData[exp][scan]['attrs']['paramtree/Scan/Step Size']*frames
    image = np.zeros((frames, pixels))
    for f in range(frames):
        image[f] = (0.45/0.5)*procData[exp][scan]['RawSpecList'][f]/rawData[exp][scan]['attrs']['paramtree/Spectrometer Camera/Exposure']
    plt.figure()
    plt.grid(b=None)
    plt.xlabel('Scan distance ($\mu$m) $\longrightarrow$')
    plt.ylabel('Pixels (spectrum axis) $\longrightarrow$')
    plt.title(exp + '/' + scan)
    plt.imshow(np.flip(np.rot90(image, 1),0), cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                interpolation='nearest', extent=(0, scanRange, pixels, 0))
    plt.colorbar(label='Photons/s', shrink=0.5)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.set_aspect(scanRange/pixels/2)
    # ax.axes.get_xaxis().set_visible(False)

    if saveFig:
        plt.savefig(path + exp + '/' + scan + '/' \
            + exp + '_' + scan + r'_spectrograph.png', \
            transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot calibration spectrograph (binned)
def plotCalSpectrograph(path, rawData, procData, exp, scan, saveFig=False):
    calFreq = rawData[exp][scan]['CalFreq']
    calFrames = calFreq.shape[0]
    pixels = procData[exp][scan]['CalSpecList'].shape[1]
    image = np.zeros((calFrames, pixels))
    for f in range(calFrames):
        image[f] = (0.45/0.5)*procData[exp][scan]['CalSpecList'][f]/rawData[exp][scan]['attrs']['paramtree/Spectrometer Camera/Ref. Exposure']
    plt.figure()
    plt.grid(b=None)
    plt.xlabel('Calibration frequency set-point (GHz) $\longrightarrow$')
    plt.ylabel('Pixels (spectrum axis) $\longrightarrow$')
    plt.title(exp + '/' + scan + '\n' + 'EOM range: %.3f to %.3f GHz' \
        %(rawData[exp][scan]['CalFreq'][0], rawData[exp][scan]['CalFreq'][-1]))
    plt.imshow(np.flip(np.rot90(image, 1),0), cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                interpolation='nearest') # extent=(calFreq[0], calFreq[-1], pixels, 0)
    plt.xticks(ticks=np.arange(0, calFrames), labels=np.round(calFreq, decimals=4))
    plt.colorbar(label='Photons/s', shrink=0.5)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.set_aspect(calFrames/pixels/2)
    # ax.set_aspect((calFreq[-1]-calFreq[0])/pixels/2)
    # ax.axes.get_xaxis().set_visible(False)

    if saveFig:
        plt.savefig(path + exp + '/' + scan + '/' \
            + exp + '_' + scan + r'_calSpectrograph.png', \
            transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot Brillouin spectra (data)
def plotSpectra(path, rawData, procData, exp, scan, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'Spectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'Spectra' + '/')
    frames = rawData[exp][scan]['attrs']['paramtree/Scan/Frame Number']
    for i in range(frames):
        plt.figure()
        plt.plot(procData[exp][scan]['FittedSpect'][i], '-', linewidth=2)
        plt.plot(procData[exp][scan]['RawSpecList'][i], 'ko', markersize=4)
        plt.title(exp + '/' + scan + ': frame %d \n Fitted BFS = %.3f ± %.3f GHz \n Fitted LW = %.3f ± %.3f GHz'\
                  %(i, procData[exp][scan]['FreqList'][i], \
                  procData[exp][scan]['FreqList_sig'][i], \
                  procData[exp][scan]['LinewidthList'][i], \
                  procData[exp][scan]['LinewidthList_sig'][i]))
        plt.xlabel('Pixels')
        plt.ylabel('Counts')
        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/'  + 'Spectra' + '/' \
                + exp + '_' + scan + r'_frame_%d.png' %i, \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot calibration spectra
def plotCalSpectra(path, rawData, procData, exp, scan, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'CalSpectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'CalSpectra' + '/')
    calFrames = rawData[exp][scan]['CalFreq'].shape[0]
    for l in range(calFrames):
        plt.figure()
        plt.plot(procData[exp][scan]['FittedCalSpect'][l], '-', linewidth=2)
        plt.plot(procData[exp][scan]['CalSpecList'][l], 'ko', markersize=4)
        plt.title(exp + '/' + scan + ': calFrame %d \n EOM setting: %.3f GHz' \
            %(l, rawData[exp][scan]['CalFreq'][l]))
        plt.xlabel('Pixels')
        plt.ylabel('Counts')
        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'CalSpectra' + '/' \
                + exp + '_' + scan + r'_calFrame_%d.png' %l, \
                transparent=True, bbox_inches='tight', pad_inches=0.01)
                

#### Plot A-line (depth) scan
def plotAline(path, rawData, procData, exp, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + 'A-lines' + '/'):
            os.makedirs(path + exp + '/' + 'A-lines' + '/')
    
    # Scans to plot
    scans = list(rawData[exp].keys())
    if scans[0]=='name':
        scans = scans[1:]

    for i, scan in enumerate(scans):
        fig, ax = plt.subplots(nrows=2, ncols=1, tight_layout=True)
        fig.set_size_inches(4, 5, True)
        
        steps = rawData[exp][scan]['attrs']['paramtree/Scan/Step Size']*np.arange(procData[exp][scan]['FreqList'].shape[0]).astype(float)
        signal = np.amax(procData[exp][scan]['FittedSpect'], 1)
        nonanIdx = ~np.isnan(procData[exp][scan]['FreqList'])
        aline = procData[exp][scan]['FreqList'][nonanIdx]  # Remove NaNs
        linew = procData[exp][scan]['LinewidthList'][nonanIdx]
        r_squared = procData[exp][scan]['FittedSpect_R^2'][nonanIdx]  # Remove NaNs
        steps = steps[nonanIdx]  # Remove NaNs
        signal = signal[nonanIdx]  # Remove NaNs

        # Eliminate points in air using signal threshold
        sigThresh = 0.2*np.mean(signal[np.argsort(signal)[-2]]) #0.5
        r_sqThresh = 0.9 #0.99
        startIdx = np.where( (signal>sigThresh) & (r_squared>r_sqThresh) )
        aline_filtered = aline[startIdx[0][0]:]
        linew_filtered = linew[startIdx[0][0]:]   
        steps_filtered = steps[startIdx[0][0]:]
        signal_filtered = signal[startIdx[0][0]:]
        r_squared_filtered = r_squared[startIdx[0][0]:]
        # Another filter to remove deeper points with too little signal
        sigIdx = np.where( (signal_filtered>0.01*sigThresh) & (r_squared_filtered>0.9) )
        
        # Plot Brillouin frequency shift vs. depth
        ax[0].plot(steps, aline, 'k.', markersize=8)
        ax[0].plot(steps_filtered[sigIdx], aline_filtered[sigIdx], \
                   'r.', markersize=8)
        # ax[0].set_xlim(0, 1000)
        ax[0].set_ylim(5.6, 6.4)
        ax[0].set_xlabel('Depth ($\mu$m)')
        ax[0].set_ylabel('Brillouin frequency shift (GHz)')
        ax[0].set_title('Scan %d/%d' %(i+1, len(scans)))
        
        # Plot Brillouin peak width vs. depth
        ax[1].plot(steps, linew, 'k.', markersize=8)
        ax[1].plot(steps_filtered[sigIdx], linew_filtered[sigIdx], \
                   'b.', markersize=8)
        # ax[1].set_xlim(0, 1000)
        ax[1].set_ylim(0.4, 1.2)
        ax[1].set_xlabel('Depth ($\mu$m)')
        ax[1].set_ylabel('Brillouin peak width (GHz)')
        
        if saveFig:
            plt.savefig(path + exp + '/' + 'A-lines' + '/' \
                + exp + '_' + scan + r'.png', \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot XZ frequency shift map from A-lines
def plotXZmap(path, rawData, procData, exp, saveFig=False):
    # Import standard Brillouin colormap
    # colors = np.genfromtxt('colormap.txt', delimiter='\t')
    # colormap = LinearSegmentedColormap.from_list('brillouin', colors, N=200)
    colormap = 'jet'
    
    # Scans to plot
    scans = list(rawData[exp].keys())
    if scans[0]=='name':
        scans = scans[1:]
    
    zOriginPos = 12160
    stepSize = rawData[exp]['Scan_100']['attrs']['paramtree/Scan/Step Size']
    xrange = 20*len(scans)
    yrange = stepSize*(15 + rawData[exp]['Scan_100']['attrs']['paramtree/Scan/Frame Number'])
    alineArr = np.full([len(scans), 15 + rawData[exp]['Scan_100']['attrs']['paramtree/Scan/Frame Number']], np.nan)

    for i, scan in enumerate(scans):
        alineRaw = np.copy(procData[exp][scan]['FreqList'])
        signal = np.amax(procData[exp][scan]['FittedSpect'], 1)
        signal[np.isnan(signal)] = 0
        sigThresh = 0.93*np.mean(signal[np.argsort(signal)[-2]]) #0.5
        startIdx = np.where(signal>sigThresh)[0][0]
        alineRaw[:startIdx] = np.nan
        zStartPos = rawData[exp][scan]['attrs']['paramtree/Motor/Current Location']
        zStartOffset = round((zStartPos - zOriginPos)/rawData[exp][scan]['attrs']['paramtree/Scan/Step Size'])
        r_squared = procData[exp][scan]['FittedSpect_R^2']
        alineRaw[r_squared<0.76] = np.nan
        alineArr[i, zStartOffset:zStartOffset+len(alineRaw)] = alineRaw
    
    ### Plot frequency shift X-Z map
    plt.matshow(np.rot90(alineArr, k=3), cmap=colormap, \
                vmin=5.17, vmax=6.08, origin='lower', interpolation=None, \
                extent=(0, xrange, -120, -120+yrange))
    # plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
    plt.ylim(0, 1100)
    plt.xlabel(r'x (${\rm \mu m}$) $\longrightarrow$')
    plt.ylabel(r'z (${\rm \mu m}$) $\longrightarrow$')
    plt.colorbar(label='Brillouin frequency shift [GHz]', shrink=0.4)
    plt.grid(b=None)
    ax = plt.gca()
    # ax.axis('off')
    ax.autoscale_view
    ax.xaxis.set_ticks_position('bottom')

    if saveFig:
        plt.savefig(path + exp + '/' \
            + exp + '_' + scan + r'_xz_map.png', \
            transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot monitor camera images for single scan
def plotBrightfield(path, rawData, exp, scan, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'Monitor images' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'Monitor images' + '/')
    frames = rawData[exp][scan]['attrs']['paramtree/Scan/Frame Number']
    for i in range(frames):
        image = np.rot90(rawData[exp][scan]['CMOSImage'][i], -1, (1,0))
        plt.figure()
        plt.grid(b=None)
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest')
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/'  + 'Monitor images' + '/' \
                + exp + '_' + scan + r'_frame_%d.png' %i, \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


