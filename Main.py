#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 5 19:07:03 2021

@author: ame38
"""
import DataImport
import Plotting
import os

path = '/Users/ame38/Downloads/2022-05-26/'
session = 'Porcine_globe'

fullPath = path + session + '/'
if not os.path.exists(fullPath):
    os.makedirs(fullPath)

# Import data and fit spectra:
filename = path + session + '.hdf5'
(rawData, procData) = DataImport.importData(filename, verbose=False)

# Create folders to save plots:
for exp in list(rawData.keys()):
    if not os.path.exists(fullPath + exp + '/'):
        os.makedirs(fullPath + exp + '/')

#%%
###### Plot A-line (depth) scans for a single exp ######
exp = 'Exp_0'

# Save figures?
saveFig = True

Plotting.plotAline(fullPath, rawData, procData, exp, saveFig)

#%%
###### Plot XZ frequency shift map from A-lines ######
exp = 'Exp_1'

# Save figures?
saveFig = False

Plotting.plotXZmap(fullPath, rawData, procData, exp, saveFig)

#%%
###### Plot Brillouin spectra (data) for a single scan ######
exp = 'Exp_1'
scan = 'Scan_50'

# Save figures?
saveFig = True
if saveFig:    
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotSpectra(fullPath, rawData, procData, exp, scan, saveFig)

#%%
###### Plot raw spectrum images (from Andor camera) for a single scan ######
exp = 'Exp_10'
scan = 'Scan_2'

# Save figures?
saveFig = True
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotAndorImage(fullPath, rawData, exp, scan, saveFig)

#%%
###### Plot spectrograph (binned) for a single scan ######
exp = 'Exp_0'
scan = 'Scan_0'

# Save figures?
saveFig = True
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotSpectrograph(fullPath, rawData, procData, exp, scan, saveFig)

#%%
###### Plot calibration spectra for a single scan ######
exp = 'Exp_0'
scan = 'Scan_8'

# Save figures?
saveFig = False
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotCalSpectra(fullPath, rawData, procData, exp, scan, saveFig)

#%%
###### Plot calibration spectrum images (from Andor camera) for a single scan ######
exp = 'Exp_0'
scan = 'Scan_8'

# Save figures?
saveFig = True
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotCalAndorImage(fullPath, rawData, exp, scan, saveFig)

#%%
###### Plot calibration spectrograph (binned) for a single scan ######
exp = 'Exp_0'
scan = 'Scan_10'

# Save figures?
saveFig = True
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotCalSpectrograph(fullPath, rawData, procData, exp, scan, saveFig)

#%%
###### Plot monitor camera images for single scan ######
exp = 'Exp_1'
scan = 'Scan_75'

# Save figures?
saveFig = False
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotBrightfield(path, rawData, exp, scan, saveFig)
