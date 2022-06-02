#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 6 15:34:54 2021

@author: ame38
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter
from lmfit import Model
from lmfit import Parameters
import itertools

#### Fit Brillouin spectrum,
# sline is the data (counts) for the pixels on the spectral line,
# ftol and xtol are fitting tolerances (adjust for speed vs. accuracy)
# rSq_thresh is a minimum r^2 threshold (to assess fit quality)
def fitSpectrum(sline, win_fac=0, rSq_thresh = 0.5, xtol=1e-6, ftol=1e-6, verbose=False):
    # Find peak locations:
    slineMax = np.amax(np.abs(sline))
    prominence = 0.3*slineMax # 0.05*np.amax(posSline)
    wlen = 5*prominence
    pk_ind, pk_info = find_peaks(sline, prominence=prominence, width=2, \
        height=50, rel_height=0.5, wlen=wlen)
    pk_wids = 0.5*pk_info['widths']
    pk_hts = pk_info['peak_heights']
    
    # Remove non-sensical peaks:
    pk_ind = pk_ind[pk_hts>0]
    pk_wids = pk_wids[pk_hts>0]
    pk_hts = pk_hts[pk_hts>0]
    pk_ind = pk_ind[pk_hts<2*slineMax]
    pk_wids = pk_wids[pk_hts<2*slineMax]
    pk_hts = pk_hts[pk_hts<2*slineMax]
    pk_ind = pk_ind[pk_wids>0]
    pk_wids = pk_wids[pk_wids>0]
    pk_hts = pk_hts[pk_wids>0]
    pk_ind = pk_ind[pk_wids<50]
    pk_wids = pk_wids[pk_wids<50]
    pk_hts = pk_hts[pk_wids<50]

    # Check for no peaks:
    if len(pk_ind)<1:
        if verbose:
            print('[DataFitting] Warning: Too few peaks in spectrum')
        interPeaksteps = np.nan
        linewidth = np.nan
        interPeaksteps_sig = np.nan
        linewidth_sig = np.nan
        fittedSpect = np.nan*np.ones(sline.shape)
        fittedSpect_rSq = np.nan
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect, fittedSpect_rSq)

    # Check for single peak (i.e. overlapping peaks):
    if len(pk_ind)==1:
        if verbose:
            print('[DataFitting] Warning: Potentially overlapping peaks')
        pk_ind = np.array([pk_ind[0]-2, pk_ind[0]+2])
        pk_wids = np.array([pk_wids[0]-2, pk_wids[0]+2])
        pk_hts = np.array([pk_hts[0]-2, pk_hts[0]+2])

    # Check for more than 2 peaks:
    if len(pk_ind)>2:
        if verbose:
            print('[DataFitting] Warning: Potentially > 2 peaks')

        # Remove significantly smaller peaks
        pk_srt = np.argsort(pk_hts)
        pk_select = pk_hts/np.nanmean(pk_hts[pk_srt[-2:]])>0.2
        pk_ind = pk_ind[pk_select]
        pk_wids = pk_wids[pk_select]
        pk_hts = pk_hts[pk_select]
        
        # If 3 or 5 peaks remain, remove asymmetric peak
        if len(pk_ind)==3 or len(pk_ind)==5:
            pk_ctr_dist = pk_ind - np.floor(0.5*sline.shape[0])
            pairs = np.array([(a,b) for a,b in itertools.combinations(np.arange(0, len(pk_ind)), 2)])
            sums = np.array([pk_ctr_dist[a]+pk_ctr_dist[b] for (a,b) in pairs])
            sum_srt = np.argsort(np.abs(sums))
            if len(pk_ind)==3:
                pk_select = np.unique(pairs[sum_srt[:1]])
            else:
                pk_select = np.unique(pairs[sum_srt[:2]])
                while len(pk_select)<4:
                    sum_srt = np.delete(sum_srt, 1)
                    pk_select = np.unique(pairs[sum_srt[:2]])
            pk_ind = pk_ind[pk_select]
            pk_wids = pk_wids[pk_select]
            pk_hts = pk_hts[pk_select]

    if len(pk_ind)==4:
        if verbose:
            print('[DataFitting] Warning: Two pairs of peaks detected')
        ### Fit spectrum to 4-Lorentzian model:
        # Starting guesses for fit:
        p0 = [pk_hts[0], pk_ind[0], pk_wids[0], \
              pk_hts[1], pk_ind[1], pk_wids[1], \
              pk_hts[2], pk_ind[2], pk_wids[2], \
              pk_hts[3], pk_ind[3], pk_wids[3], \
              np.amin(sline)]
        
        pix = np.arange(0, sline.shape[0]) # Pixel number
        
        # Create boolean mask to filter out points far from the peaks:
        pk_mask = np.array(0*sline, dtype=bool)
        if win_fac > 0:
            pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
            pk_mask[(pk_ind[1] - win_fac*pk_wids[1]).astype(int):(pk_ind[1] + win_fac*pk_wids[1]).astype(int)]=True
            pk_mask[(pk_ind[2] - win_fac*pk_wids[2]).astype(int):(pk_ind[2] + win_fac*pk_wids[2]).astype(int)]=True
            pk_mask[(pk_ind[3] - win_fac*pk_wids[3]).astype(int):(pk_ind[3] + win_fac*pk_wids[3]).astype(int)]=True
        else:
            pk_mask[:]=True # Do not use mask
        
        # Fit spectrum to 4-Lorentzian model
        try:
            popt, pcov = curve_fit(_4Lorentzian, pix[pk_mask], sline[pk_mask], \
                              p0=p0, ftol=ftol, xtol=xtol)
            psig = np.sqrt(np.diag(pcov))
            
            # Remove any non-sensical peaks
            pk_ht_ind = np.array([0, 3, 6, 9])
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] > 0]
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] < 2*slineMax]
            pk_ht_ind = pk_ht_ind[np.abs(popt[pk_ht_ind+2]) < 50]
            
            # Find larger pair of peaks
            pk_ht_srt = np.argsort(popt[pk_ht_ind])
            pk_1_ind = pk_ht_ind[pk_ht_srt[-1]]
            pk_2_ind = pk_ht_ind[pk_ht_srt[-2]]
            interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
            linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
            interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
            linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
            fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
                                        popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[12])
            # fittedSpect = _4Lorentzian(pix[pk_mask], \
            #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
            #                 popt[5], popt[6], popt[7], popt[8], popt[9], \
            #                 popt[10], popt[11], popt[12])
            residuals = sline - fittedSpect
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sline-np.mean(sline))**2)
            fittedSpect_rSq = 1 - (ss_res / ss_tot)
            
        except:
            if verbose:
                print('[DataFitting] Warning: Fitting 4-peak spectrum failed')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            fittedSpect = np.nan*np.ones(sline.shape)
            fittedSpect_rSq = np.nan
            return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect, fittedSpect_rSq)
        
        # R^2 quality check
        if fittedSpect_rSq < rSq_thresh:
            if verbose:
                print('[DataFitting] Warning: Fitting 4-peak spectrum failed - low R^2')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            fittedSpect = np.nan*np.ones(sline.shape)
            fittedSpect_rSq = np.nan
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect, fittedSpect_rSq)
                
    ### Fit spectrum to 2-Lorentzian model:
    # Starting guesses for fit:
    p0 = [pk_hts[0], pk_ind[0], pk_wids[0], \
          pk_hts[1], pk_ind[1], pk_wids[1], \
          np.amin(sline)]
    
    pix = np.arange(0, sline.shape[0]) # Pixel number
    
    # Create boolean mask to filter out points far from the peaks:
    pk_mask = np.array(0*sline, dtype=bool)
    if win_fac > 0:
        pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
        pk_mask[(pk_ind[1] - win_fac*pk_wids[1]).astype(int):(pk_ind[1] + win_fac*pk_wids[1]).astype(int)]=True
    else:
        pk_mask[:]=True # Do not use mask

    # Fit spectrum to 2-Lorentzian model:
    try:
        popt, pcov = curve_fit(_2Lorentzian, pix[pk_mask], sline[pk_mask], \
                          p0=p0, ftol=ftol, xtol=xtol)
        psig = np.sqrt(np.diag(pcov))
        
        # Check for peak height imbalance (indicating overlapping peaks)
        if popt[0]>3*popt[3]:
            if verbose:
                print('[DataFitting] Warning: Overlapping peaks detected')
            interPeaksteps = 0
            linewidth = np.abs(popt[2])
            interPeaksteps_sig = psig[1]
            linewidth_sig = psig[2]
            fittedSpect = _2Lorentzian(pix, popt[0], popt[1], popt[2], 0, 1, 1, popt[6])
            residuals = sline - fittedSpect
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sline-np.mean(sline))**2)
            fittedSpect_rSq = 1 - (ss_res / ss_tot)

        elif popt[3]>3*popt[0]:
            if verbose:
                print('[DataFitting] Warning: Overlapping peaks detected')
            interPeaksteps = 0
            linewidth = np.abs(popt[5])
            interPeaksteps_sig = psig[4]
            linewidth_sig = psig[5]
            fittedSpect = _2Lorentzian(pix, 0, 1, 1, popt[3], popt[4], popt[5], popt[6])
            residuals = sline - fittedSpect
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sline-np.mean(sline))**2)
            fittedSpect_rSq = 1 - (ss_res / ss_tot)
            
        else:
            interPeaksteps = np.abs(popt[4] - popt[1])
            linewidth = 0.5*(np.abs(popt[2]) + np.abs(popt[5])) # Mean linewidth of 2 peaks
            interPeaksteps_sig = np.sqrt(psig[4]**2 + psig[1]**2)
            linewidth_sig = 0.5*np.sqrt(psig[2]**2 + psig[5]**2)
            fittedSpect = _2Lorentzian(pix, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
            residuals = sline - fittedSpect
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sline-np.mean(sline))**2)
            fittedSpect_rSq = 1 - (ss_res / ss_tot)
            
    except:
        if verbose:
            print('[DataFitting] Warning: Fitting 2-peak spectrum failed')
        interPeaksteps = np.nan
        linewidth = np.nan
        interPeaksteps_sig = np.nan
        linewidth_sig = np.nan
        fittedSpect = np.nan*np.ones(sline.shape)
        fittedSpect_rSq = np.nan
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect, fittedSpect_rSq)
    
    # R^2 quality check
    if fittedSpect_rSq < rSq_thresh:
        if verbose:
            print('[DataFitting] Warning: Fitting 2-peak spectrum failed - low R^2')
        interPeaksteps = np.nan
        linewidth = np.nan
        interPeaksteps_sig = np.nan
        linewidth_sig = np.nan
        fittedSpect = np.nan*np.ones(sline.shape)
        fittedSpect_rSq = np.nan
    return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect, fittedSpect_rSq)


#### Fit calibration curve to determine SD and FSR
def fitCalCurve(pxDist, freq, xtol=1e-6, ftol=1e-6, verbose=False):
    # Starting guesses for fit:
    p0 = [0.127, 21.5]
    # Remove NaNs before fitting
    freq = freq[~np.isnan(pxDist)]
    pxDist = pxDist[~np.isnan(pxDist)]
    try:
        popt, pcov = curve_fit(_Linear, pxDist, freq, p0=p0, ftol=ftol, xtol=xtol)
    except:
        if verbose:
            print('[DataFitting] Warning: Fitting calibration curve failed')
        SD = np.nan
        FSR = np.nan
        SD_sig = np.nan
        FSR_sig = np.nan
        return (SD, FSR, SD_sig, FSR_sig)
    residuals = freq - _Linear(pxDist, popt[0], popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((freq-np.mean(freq))**2)
    r_squared = 1 - (ss_res / ss_tot)
    # print('[DataFitting] Calibration curve fitting r^2 = %.4f' %r_squared)
    if r_squared < 0.9:
        if verbose:
            print('[DataFitting] Warning: Fitting calibration curve failed - low R^2')
        SD = np.nan
        FSR = np.nan
        SD_sig = np.nan
        FSR_sig = np.nan
        return (SD, FSR, SD_sig, FSR_sig)
    else:
        psig = np.sqrt(np.diag(pcov))
        SD = popt[0]
        FSR = popt[1]
        SD_sig = psig[0]
        FSR_sig = psig[1]
    return (SD, FSR, SD_sig, FSR_sig)


#### Segment Brillouin A-line (depth scan) into air / cornea / aq humor
def fitAline(aline, r_squared, signal, xtol=1e-6, ftol=1e-6, maxfev=1000):
    steps = np.arange(aline.shape[0]).astype(float)
    nonanIdx = ~np.isnan(aline)
    aline = aline[nonanIdx]  # Remove NaNs
    steps = steps[nonanIdx]  # Remove NaNs
    signal = signal[nonanIdx]  # Remove NaNs
    r_squared = r_squared[nonanIdx]  # Remove NaNs
    
    # Eliminate points in air using signal threshold
    sigThresh = 0.3*np.mean(signal[np.argsort(signal)[-4:-1]]) #0.5, 0.3
    startIdx = np.where( (signal>sigThresh) & (r_squared>0.85) ) #0.99, 0.85
    aline = aline[startIdx[0][0]:]
    steps = steps[startIdx[0][0]:]
    signal = signal[startIdx[0][0]:]
    
    # Remove points with freq values out of range
    cleanIdx = (aline > 5) & (aline < 6)
    
    # Check there are enough points in A-line
    if np.sum(cleanIdx==True) < 3:
        print('[DataFitting] Too few valid data points. Could not fit A-line')
        stromaBFS = np.nan
        stromaIdx = []
        popt = [0, 0, 0, 0]
        return (stromaBFS, stromaIdx, popt)
    
    # Check if scan only includes points in cornea (using std dev)
    if np.std(aline[cleanIdx]) < 0.05:
        print('[DataFitting] A-line likely contains only points in stroma')
        stromaBFS = np.mean(aline[cleanIdx])
        stromaIdx = steps[cleanIdx].astype(int)
        popt = [0, 0, 0, 0]
        return (stromaBFS, stromaIdx, popt)
    
    aline = aline[cleanIdx]    
    steps = steps[cleanIdx]
    signal = signal[cleanIdx]
    
    # Remove 1st/2nd points if outlier
    while np.abs(aline[1]-aline[0]) > 0.04:
        aline = aline[1:]
        steps = steps[1:]
        signal = signal[1:]
        if aline.shape[0]<3:
            break
    while np.abs(aline[2]-aline[0]) > 0.04:
        aline = aline[2:]
        steps = steps[2:]
        signal = signal[2:]
        if aline.shape[0]<3:
            break
    
    # Smooth with Gaussian filter, then take derivative
    blurred = gaussian_filter(aline, sigma=1)
    deriv = np.gradient(blurred, steps)
    edgeIdx = np.argmin(deriv)
    sortAline = np.argsort(blurred)
    if len(blurred) > 8:
        stromaGuess = np.mean(blurred[sortAline[-4:]])
        aqGuess = np.mean(blurred[sortAline[:4]])
        # print('[DataFitting] stromaGuess = %.2f GHz' %stromaGuess)
        # print('[DataFitting] aqGuess = %.2f GHz' %aqGuess)
    else:
        stromaGuess = np.amax(blurred, axis=0)
        aqGuess = np.amin(blurred, axis=0)
    # Estimate width of transition region from stroma to aqueous humor
    widthGuess = (aqGuess - stromaGuess)/deriv[edgeIdx]
    endGuess = steps[edgeIdx]-0.5*widthGuess
    # Starting guesses for fit:
    p0 = [endGuess, stromaGuess, widthGuess, aqGuess]
    
    # Fit A-line shape to piecewise function
    try:
        popt, pcov = curve_fit(_AlineShape, steps, aline, p0=p0, ftol=ftol, xtol=xtol)
        residuals = aline - _AlineShape(steps, popt[0], popt[1], popt[2], popt[3])
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((aline-np.mean(aline))**2)
        r_squared = 1 - (ss_res / ss_tot)
        # psig = np.sqrt(np.diag(pcov))
        idxEnd = (np.abs(steps - popt[0])).argmin()
        strMask = (steps < steps[idxEnd]) & (abs(aline - stromaGuess) < 0.1)
        strIdx = np.where(strMask)[0]
    except:
        print('[DataFitting] Could not fit A-line')
        stromaBFS = np.nan
        stromaIdx = []
        popt = [0, 0, 0, 0]
        return (stromaBFS, stromaIdx, popt)
    
    # # Only include anterior 1/2 of stroma
    # endIdx = int(np.round(0.5*np.sum(strIdx==True)))
    # strIdx[endIdx:] = False
    
    # Isolate stroma points
    if len(strIdx) > 0:
        strPts = aline[strIdx]
        # Remove 1st/2nd point if outlier:
        while np.abs(strPts[1]-strPts[0]) > 0.04:
            strIdx = strIdx[1:]
            strPts = aline[strIdx]
            if strPts.shape[0]<3:
                break
        while np.abs(strPts[2]-strPts[0]) > 0.04:
            strIdx = strIdx[2:]
            strPts = aline[strIdx]
            if strPts.shape[0]<3:
                break
    else:
        print('[DataFitting] Too few stroma points')
        stromaBFS = np.nan
        stromaIdx = []
        return (stromaBFS, stromaIdx, popt)
    
    # Calculate mean BFS in segmented stroma pts
    try:
        stromaBFS = np.nanmean(strPts)
        stromaIdx = steps[strIdx].astype(int)
    except:
        print('[DataFitting] Could not fit A-line')
        stromaBFS = np.nan
        stromaIdx = []
        return (stromaBFS, stromaIdx, popt)
    return (stromaBFS, stromaIdx, popt)


#### Fit Brillouin angle-dependence (anisotropy) curve w/ 2 variables
def fitAngleDep2Var(angle, freq):
    pars = Parameters()
    pars.add('delta', value=0.0, min=0, max=1/3) #min=0.0, max=0.2
    pars.add('epsilon', value=0.1, min=0.0, max=0.5)
    if freq[angle==0].shape[0] < 1:
        freq0 = np.mean(freq[np.argsort(angle)[:3]])
    else:
        freq0 = np.nanmean(freq[angle==0])
    try:
        model2 = Model(_Anisotropy2)
        result = model2.fit(freq/freq0, pars, x=angle, method='nelder')
        print(result.fit_report())
    except:
        print('[DataFitting/fitAngleDep2Var] Could not fit angle dependence')
        result = None
    return result


#### Fit Brillouin angle-dependence (anisotropy) curve w/ 3 variables
def fitAngleDep3Var(angle, freq):
    pars = Parameters()
    pars.add('delta', value=0.0, min=0, max=1/3) #min=0.0, max=0.2
    pars.add('epsilon', value=0.1, min=0.0, max=0.5)
    pars.add('omega0', value=5.5, min=5.0, max=6.0)
    try:
        model3 = Model(_Anisotropy3)
        result = model3.fit(freq, pars, x=angle, method='nelder')
        print(result.fit_report())
    except:
        print('[DataFitting/fitAngleDep3Var] Could not fit angle dependence')
        result = None
    return result


############  Fitting functions  ############ 

def _2Lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + (amp2*wid2**2/((x-cen2)**2+wid2**2)) \
            + offs

def _4Lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + (amp2*wid2**2/((x-cen2)**2+wid2**2)) \
            + (amp3*wid3**2/((x-cen3)**2+wid3**2)) \
            + (amp4*wid4**2/((x-cen4)**2+wid4**2)) \
            + offs

def _Linear(x, sd, fsr):
    return 0.5*fsr - 0.5*sd*x

def _AlineShape(x, endStroma, stromaBFS, width, aqBFS):
    condlist = [ x < endStroma,
                (x >= endStroma) & (x < endStroma + width),
                x >= endStroma + width
                ]
    funclist = [lambda x: stromaBFS,
                lambda x: x*(aqBFS-stromaBFS)/width + stromaBFS \
                    - endStroma*(aqBFS-stromaBFS)/width,
                lambda x: aqBFS
                ]
    return np.piecewise( x, condlist, funclist )

def _Anisotropy2(x, delta, epsilon):
    return np.sqrt(1 + \
                2*delta*epsilon*np.sin(np.deg2rad(x))**2*np.cos(np.deg2rad(x))**2 \
                + epsilon*np.sin(np.deg2rad(x))**4)

def _Anisotropy3(x, delta, epsilon, omega0):
    return omega0*np.sqrt(1 + \
                2*delta*epsilon*np.sin(np.deg2rad(x))**2*np.cos(np.deg2rad(x))**2 \
                + epsilon*np.sin(np.deg2rad(x))**4)
