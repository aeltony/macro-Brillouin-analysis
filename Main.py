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
# path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Data/' + \
#         'Sensitivity/2022-02-11/'
# path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Data/' + \
#         'Porcine limbus and sclera/2022-02-10/'
# session = 'Cuvette_sensitivity'
# session = 'Exposure_time_water'
# session = 'Calibration_example'
# session = 'Rb_vapor_extinction'
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

#%%
###### Fit A-line shapes and plot angle-dep.for porcine cornea ######
import CornealAnisotropy

# Exp's to plot
exps = list(rawData.keys())
exps = exps[1:]

CornealAnisotropy.fitPorcineAngle(fullPath, rawData, procData, exps)

#%%
###### Plot aggregate porcine cornea angle data ######
import CornealAnisotropy

fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/' + \
            'Projects/Anisotropy of cornea/2021 New Data/Porcine angle/' + \
            'Porcine_angle_data_3var.txt'

CornealAnisotropy.plotAggregateAngleData(fullPath)

#%%
###### Fit A-line shapes for all scans in Brillouin map (single exp) ######
import CornealAnisotropy

# Exp to plot
exp = 'Exp_1'

CornealAnisotropy.fitPorcineMap(fullPath, rawData, procData, exp)

#%%
###### Plot aggregate porcine globe map data ######
import CornealAnisotropy

fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/' + \
        'Projects/Anisotropy of cornea/2021 New Data/Porcine maps/' + \
        'Porcine_map_data_3var.txt'

CornealAnisotropy.plotAggregateMapData(fullPath)

#%%
###### Plot human map data ######
import CornealAnisotropy

path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' + \
        'Anisotropy of cornea/2021 New Data/Human maps/'

# session = 'DB-right' # Figure 4A
# session = 'P10002OS-Scaled' # Figure 4B
# session = 'P10201-OS' # Figure 4C
session = 'P10200-OS' # Figure 4D

CornealAnisotropy.fitHumanMap(path, session)

#%%
###### Calculate aggregate human map data ######
import CornealAnisotropy

fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' + \
        'Anisotropy of cornea/2021 New Data/Human maps/' + \
        'Human_map_data_3var.txt'

CornealAnisotropy.plotAggregateHumanData(fullPath)

