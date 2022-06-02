#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 7 17:15:35 2021

@author: ame38
"""
import h5py
import numpy as np
import DataFitting
import os

#### Import all data from hdf5 file (output by DAQ software)
# rawData is dictionary of raw data and parameters used
# procData is dictionary of processed data (i.e. fitted + calibrated spectra)
def importData(filename, verbose=False):
    # Create text file for saving exp notes
    name = os.path.splitext(os.path.basename(filename))[0]
    path = os.path.dirname(filename) + '/'
    f_notes = open(path + name + '/' + name + '.txt', 'w')
    f_notes.write('Data Summary: \n')
    # Open data file and load hdf5 structure to dictionary
    f = h5py.File(filename, 'r')
    rawData = {}
    procData = {}
    win_fac = 0 # Number of points either side of peak to include in fit, 0 = no window
    for e in list(f.keys()):
        rawData[e] = {}
        procData[e] = {}
        if f[e].attrs.__contains__('name'):
            rawData[e]['name'] = f[e].attrs['name']
            print(e + ': ' + f[e].attrs['name'])
            f_notes.write('\n' + e + ': ' + f[e].attrs['name'] + '\n')
        scans = list(f[e].keys())
        scanNums = [s[5:] for s in scans]
        scanNums.sort(key = int)
        for n in scanNums:
            s = 'Scan_' + str(n)
            rawData[e][s] = {}
            procData[e][s] = {}
            rawData[e][s]['attrs'] = dict(f[e][s].attrs.items())
            for k in list(f[e][s].keys()):
                rawData[e][s][k] = np.array(f[e][s][k])
            if rawData[e][s]['attrs']['note']:
                print('Processing ' + e + '/' + s + ', ' + rawData[e][s]['attrs']['note'])
                f_notes.write(rawData[e][s]['attrs']['timestamp'] + ' - ' + s + ': ' + rawData[e][s]['attrs']['note'] + '\n')
            else:
                print('Processing ' + e + '/' + s)
                f_notes.write(rawData[e][s]['attrs']['timestamp'] + ' - ' + s + '\n')
            # Separate sample vs. calibration frames
            calFrames = rawData[e][s]['CalFreq'].shape[0]
            frames = rawData[e][s]['attrs']['paramtree/Scan/Frame Number']
            procData[e][s]['RawSpecList'] = np.copy(rawData[e][s]['SpecList'][:-calFrames])
            procData[e][s]['CalSpecList'] = np.copy(rawData[e][s]['SpecList'][-calFrames:])     
            # Fitting Brillouin spectra
            if procData[e][s]['RawSpecList'].shape[0] > 0:
                procData[e][s]['IntegrPhotonsList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['FreqList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['LinewidthList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['FreqList_sig'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['LinewidthList_sig'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['FittedSpect'] = np.empty(procData[e][s]['RawSpecList'].shape)
                procData[e][s]['FittedSpect_R^2'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
                procData[e][s]['FittedCalSpect'] = np.empty(procData[e][s]['CalSpecList'].shape)
                procData[e][s]['FittedCalSpect_R^2'] = np.zeros(procData[e][s]['CalSpecList'].shape[0])
                # Find SD / FSR
                procData[e][s]['PxDist'] = np.empty(rawData[e][s]['CalFreq'].shape)
                procData[e][s]['PxDist_sig'] = np.empty(rawData[e][s]['CalFreq'].shape)
                for i in range(calFrames):
                    procData[e][s]['PxDist'][i], width, \
                        procData[e][s]['PxDist_sig'][i], width_sig, \
                        procData[e][s]['FittedCalSpect'][i], \
                        procData[e][s]['FittedCalSpect_R^2'][i] = \
                        DataFitting.fitSpectrum(np.copy(procData[e][s]['CalSpecList'][i]),
                                                win_fac, 0.9, 1e-7, 1e-7, verbose)
                try:
                    procData[e][s]['SDcal'], procData[e][s]['FSRcal'], \
                        procData[e][s]['SDcal_sig'], procData[e][s]['FSRcal_sig'] = \
                        DataFitting.fitCalCurve(np.copy(procData[e][s]['PxDist']), \
                        np.copy(rawData[e][s]['CalFreq']), 1e-7, 1e-7, verbose)
                    # print('[importData] Fitted SD = %.3f ± %.5f GHz/px' \
                    #       %(procData[e][s]['SDcal'], procData[e][s]['SDcal_sig']))
                    # print('[importData] Fitted FSR = %.3f ± %.3f GHz' \
                    #       %(procData[e][s]['FSRcal'], procData[e][s]['FSRcal_sig']))
                except:
                    procData[e][s]['SDcal'] = np.nan
                    procData[e][s]['FSRcal'] = np.nan
                    if verbose:
                        print('[importData] Could not fit SD + FSR')
                # Use SD + FSR to calculate Brillouin frequency shifts
                for i in range(frames):
                    sline = np.copy(procData[e][s]['RawSpecList'][i])
                    sline = np.transpose(sline)
                    procData[e][s]['IntegrPhotonsList'][i] = \
                        (0.45/0.5)*np.sum(sline)
                    interPeakDist, width, interPeakDist_sig, width_sig, \
                        procData[e][s]['FittedSpect'][i], \
                        procData[e][s]['FittedSpect_R^2'][i] = \
                        DataFitting.fitSpectrum(sline, win_fac, 0.1, 1e-7, 1e-7, verbose)
                    procData[e][s]['FreqList'][i] = \
                        0.5*(procData[e][s]['FSRcal'] - \
                        procData[e][s]['SDcal']*interPeakDist)
                    procData[e][s]['FreqList_sig'][i] = \
                        0.5*np.sqrt( procData[e][s]['FSRcal_sig']**2 + \
                        (interPeakDist*procData[e][s]['SDcal_sig'])**2 + \
                        (interPeakDist_sig*procData[e][s]['SDcal'])**2)
                    procData[e][s]['LinewidthList'][i] = \
                        procData[e][s]['SDcal']*width
                    procData[e][s]['LinewidthList_sig'][i] = \
                        np.sqrt((procData[e][s]['SDcal']*width_sig)**2 +\
                        (procData[e][s]['SDcal_sig']*width)**2)
            else:
                if verbose:
                    print('[importData] Empty spectral line (could not fit spectrum)')
    f_notes.close()
    f.close()
    return (rawData, procData)
