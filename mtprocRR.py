# -*- coding: utf-8 -*-
"""
Created on Mon May  4 17:25:49 2020

@author: AJITHABH K. S.
Last modified: 30-05-2023

This is the main module of this package, especially written for remote reference.

It consists of functions to read time series data, detrending, calibration,
computation of fourier transform, calculation of impedance, etc.

Purpose of each functions are given as comments.
"""
import os
from scipy import signal
from scipy.fft import fft
import numpy as np
from datetime import datetime as dt
from matplotlib import pyplot as plt
from tqdm import tqdm
from scipy.signal import butter, lfilter, freqz
#
#
def readADU07e(filename):
    """
    This function is to read ADU07e raw data (.ats file)
    It will read header information and time series data from the file.

    Parameters
    ----------
    filename : It is string value that represent the path to the data file.

    Returns
    -------
    header : It is a Python dictionary having header information of the time series
        data file.
    ts : It is a Python dictionary containing arrays which are time series data for
        all five components (Ex, Ey, Hx, Hy, Hz).

    """
    with open(filename, 'rb') as f:
        header = {}
        header['length'] = np.fromfile(f, dtype=np.int16, count=1).tolist() # Header length
        header['ver'] = np.fromfile(f, dtype=np.int16, count=1).tolist() # Header version
        header['nsamples'] = np.fromfile(f, dtype=np.int32, count=1).tolist() # Number of samples
        header['sfreq'] = np.fromfile(f, dtype=np.float32, count=1).tolist() # sampling frequency, Hz
        header['start'] = np.fromfile(f,dtype=np.int32, count=1).tolist() # Start time, seconds since 1970
        header['lsb'] = np.fromfile(f,dtype=np.double, count = 1).tolist() # LSB-Value
        header['iGMTOffset'] = np.fromfile(f,dtype=np.int32, count = 1).tolist()
        header['rOrigSampleFreq'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['adu06'] = np.fromfile(f,dtype=np.int16, count = 1).tolist() # ADU serial number
        header['adu06ADB'] = np.fromfile(f,dtype=np.int16, count = 1).tolist() # ADU serial number
        header['ch_no'] = np.fromfile(f,dtype=np.int8, count = 1).tolist() # Channel number
        header['bychopper'] =  np.fromfile(f,dtype=np.int8, count = 1).tolist() # Chopper
        header['ch_type']  = np.fromfile(f,dtype=np.int8, count = 2).tolist() #channel type (Hx,Hy,...)
        header['ch_type'] = ''.join([chr(item) for item in header['ch_type']])
        header['sensor'] = np.fromfile(f,dtype=np.int8, count = 6).tolist()
        header['sensor'] = ''.join([chr(item) for item in header['sensor']])
        header['sensor_no'] = np.fromfile(f,dtype=np.int16, count = 1).tolist()  # Sensor serial number
        header['x1'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   x1 coordinate of 1. Dipole (m)
        header['y1'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   y1 coordinate of 1. Dipole (m)
        header['z1'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   z1 coordinate of 1. Dipole (m)
        header['x2'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   x2 coordinate of 1. Dipole (m)
        header['y2'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   y2 coordinate of 1. Dipole (m)
        header['z2'] = np.fromfile(f,dtype=np.float32, count = 1).tolist()  #   z2 coordinate of 1. Dipole (m)
        header['dipol_length'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['dipol_angle'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['rProbeRes'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['rDCOffset'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['rPreGain'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['rPostGain'] =  np.fromfile(f,dtype=np.float32, count = 1).tolist()
        header['lat'] = np.fromfile(f,dtype=np.int32, count = 1).tolist()
        header['lon'] = np.fromfile(f,dtype=np.int32, count = 1).tolist()
        header['elev'] = np.fromfile(f,dtype=np.int32, count = 1).tolist()
        header['byLatLongType']  = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['byLatLongType'] = ''.join([chr(item) for item in header['byLatLongType']])
        header['byAddCoordType']  = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['byAddCoordType'] = ''.join([chr(item) for item in header['byAddCoordType']])
        header['siRefMedian'] = np.fromfile(f,dtype=np.int16, count = 1).tolist()
        header['dblNorthing'] = np.fromfile(f,dtype=np.float64, count = 1).tolist()
        header['dblEasting'] = np.fromfile(f,dtype=np.float64, count = 1).tolist()
        header['byGPSStat'] = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['byGPSStat'] = ''.join([chr(item) for item in header['byGPSStat']])
        header['byGPSAccuracy'] = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['byGPSAccuracy'] = ''.join([chr(item) for item in header['byGPSAccuracy']])
        header['iUTCOffset'] = np.fromfile(f,dtype=np.int16, count = 1).tolist()
        header['achSystemType'] = np.fromfile(f,dtype=np.int8, count = 12).tolist()
        header['achSystemType'] = ''.join([chr(item) for item in header['achSystemType']])
        header['achSurveyHeaderName'] = np.fromfile(f,dtype=np.int8, count = 12).tolist()
        header['achSurveyHeaderName'] = ''.join([chr(item) for item in header['achSurveyHeaderName']])
        header['achMeasType'] = np.fromfile(f,dtype=np.int8, count = 4).tolist()
        header['achMeasType'] = ''.join([chr(item) for item in header['achMeasType']])
        header['DCOffsetCorrValue'] = np.fromfile(f,dtype=np.float64, count = 1).tolist()
        header['DCOffsetCorrOn'] = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['InputDivOn'] = np.fromfile(f,dtype=np.int8, count = 1).tolist()
        header['bit_indicator'] = np.fromfile(f,dtype=np.int16, count = 1).tolist()
        header['achSelfTestResult'] = np.fromfile(f,dtype=np.int8, count = 2).tolist()
        header['achSelfTestResult'] = ''.join([chr(item) for item in header['achSelfTestResult']])
        header['numslice'] = np.fromfile(f,dtype=np.uint16, count = 1).tolist()
        header['siCalFreqs'] = np.fromfile(f,dtype=np.uint16, count = 1).tolist()
        f.seek(header['length'][0])
        ts = np.fromfile(f,dtype=np.int32, count = header['nsamples'][0])
        ts = ts * header['lsb'][0]
        return header, ts

def ts(path):
    """
    Function to return time series data.

    Parameters
    ----------
    path : It is a string value that is the path to the measurement folder.

    Returns
    -------
    ts : It is a Python dictionary containing arrays which are time series data for
        all five components (Ex, Ey, Hx, Hy, Hz).
    fs : It is a float value representing the sampling frequency of the measurement.
    sensor_no : It is a Python dictionary containing the sensor numbers of the 
        coils used.
    timeline : It is a list of datetime of each sample in the measurement.
    ChoppStat: It is an integer value showing chopper status (related to induction coil).
    loc : It is a Python dictionary containing latitude, longitude and elevation
        details of the selected site.

    """
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.ats' in file:
                files.append(os.path.join(r, file))
    del r,d,f,file,path
    files = sorted(files)
    [headerEx,tsEx] = readADU07e(files[0])
    [headerEy,tsEy] = readADU07e(files[1])
    [headerHx,tsHx] = readADU07e(files[2])
    [headerHy,tsHy] = readADU07e(files[3])
    [headerHz,tsHz] = readADU07e(files[4])
    dipoleNS = abs(headerEx.get('x2')[0]) + abs(headerEx.get('x1')[0]) 
    dipoleEW = abs(headerEy.get('y2')[0]) + abs(headerEy.get('y1')[0])
    ts = {}
    ts['tsEx'] = tsEx/(1*dipoleNS/1000)
    ts['tsEy'] = tsEy/(1*dipoleEW/1000)
    ts['tsHx'] = tsHx
    ts['tsHy'] = tsHy
    ts['tsHz'] = tsHz
    fs = headerHx.get('sfreq')
    sensor_no = {'Hx': headerHx.get('sensor_no')}
    sensor_no['Hy'] = headerHy.get('sensor_no')
    sensor_no['Hz'] = headerHz.get('sensor_no')
    ChoppStat = headerHz.get('bychopper')
    loc = {}
    loc['lat'] = headerEx.get('lat')[0]/1000/60/60
    loc['lon'] = headerEx.get('lon')[0]/1000/60/60
    loc['elev'] = headerEx.get('elev')[0]/100
    #Timeline
    start_time = [headerEx.get('start')[0]]
    d = dt.strptime('1970-1-1 00:00','%Y-%m-%d %H:%M')
    GPS_initial = datenum(d)
    timeline = [start_time[0]/3600/24 + GPS_initial]
    for x in range(1,np.shape(ts.get('tsEx'))[0]):
        start_time.append(start_time[x-1] + 1/fs[0])
        timeline.append(start_time[x]/3600/24 + GPS_initial)
    return ts,fs[0],sensor_no,timeline,ChoppStat[0],loc


def datenum(d):
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def FFTLength(nofsamples):
    """
    Function to select a preferred FFT length.

    Parameters
    ----------
    nofsamples : It is an integer value representing the number of samples in each measurement.

    Returns
    -------
    WindowLength : It is an integer value representing the window length.

    """
    #Based on Borah & Patro, 2015
    i = 1
    FFTs = [256]
    cFFT = 256 * (2 ** i)
    term = nofsamples / (20*i);
    FFTs.append(cFFT)
    while cFFT <= term:
        i = i + 1
        cFFT = 256 * (2 ** i)
        term = nofsamples / (20*i)
        if cFFT <= term:
            FFTs.append(cFFT)
    WindowLength = FFTs[-1]
    if WindowLength > 16384:
        WindowLength = 16384
    return WindowLength

def bandavg(ts,procinfo,tsR,procinfoR,config):
    """
    This is an important function. Here, auto and cross spectra are calculated for 
        each target frequencies and impedance values are computed for all events.

    Parameters
    ----------
    ts : It is a Python dictionary containing arrays which are time series data for
        all five components (Ex, Ey, Hx, Hy, Hz).
    procinfo : It is a Python dictionary containing few information regarding the processing.
    tsR : It is a Python dictionary containing arrays which are time series data for
        all five components (Ex, Ey, Hx, Hy, Hz) of remote site.
    procinfoR : It is a Python dictionary containing few information of remote site
        regarding the processing.
    config : It is a python dictionary containing processing parameters such as
        FFT length, Parzen window radius, Mahalanobis distance threshold values.

    Returns
    -------
    ftlist : It is an array containing the target frequency list.
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values and impedance values for all time windows at all target 
        frequencies.

    """
    WindowLength = procinfo.get('WindowLength')
    overlap = procinfo.get('overlap')
    ChoppStat = procinfo.get('ChoppStat')
    ChoppStatR = procinfo.get('ChoppStat')
    fs = procinfo.get('fs')
    notch = procinfo.get('notch')
    CalDataHx, CalDataHy, CalDataHz = getcalibrationdata(procinfo)
    CalDataRx, CalDataRy, CalDataRz = getcalibrationdata(procinfoR)
    ChoppDataHx = CalDataHx.get(str(ChoppStat))
    ChoppDataHy = CalDataHy.get(str(ChoppStat))
    ChoppDataHz = CalDataHz.get(str(ChoppStat))
    ChoppDataRx = CalDataRx.get(str(ChoppStatR))
    ChoppDataRy = CalDataRy.get(str(ChoppStatR))
    cal = {}
    # Make calibration values
    magnitude = ChoppDataHx[:,0] * ChoppDataHx[:,1]
    phase = np.radians(ChoppDataHx[:,2])
    calt_Hx = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
    magnitude = ChoppDataHy[:,0] * ChoppDataHy[:,1]
    phase = np.radians(ChoppDataHy[:,2])
    calt_Hy = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
    magnitude = ChoppDataHz[:,0] * ChoppDataHz[:,1]
    phase = np.radians(ChoppDataHz[:,2])
    calt_Hz = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
    magnitude = ChoppDataRx[:,0] * ChoppDataRx[:,1]
    phase = np.radians(ChoppDataRx[:,2])
    calt_Rx = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
    magnitude = ChoppDataRy[:,0] * ChoppDataRy[:,1]
    phase = np.radians(ChoppDataRy[:,2])
    calt_Ry = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
    del magnitude, phase
    ftlist = targetfreq(fs) # get frequency list
    ftlist = np.asarray(ftlist)
    ftlist = ftlist.reshape(-1, 1)
    dof = np.empty((np.shape(ftlist)[0],1),dtype=int)
    avgf = np.empty((np.shape(ftlist)[0],1),dtype=int)
    Ex = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Ey = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Hx = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Hy = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Hz = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Rx = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    Ry = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    ExRxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    ExRyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    EyRxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    EyRyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HxRxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HxRyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HyRxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HyRyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    ExExc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    EyEyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    ExEyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HzHxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    HzHyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    tHxHxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    tHyHyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    tHxHyc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    tHyHxc = np.empty((np.shape(ftlist)[0],procinfo.get('nstacks')),dtype=complex)
    s1 = np.empty(procinfo.get('nstacks')+1,dtype=int)
    s2 = np.empty(procinfo.get('nstacks')+1,dtype=int)
    s1[0] = 0
    s2[0] = WindowLength
    tsEx = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsEy = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsHx = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsHy = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsHz = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsRx = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    tsRy = np.empty((WindowLength,procinfo.get('nstacks')),dtype=float)
    for stack in (range(procinfo.get('nstacks'))):
        #tsEx[:,stack] = ts.get('tsEx')[s1[0]:s2[0]]
        tsEx[:,stack] = ts.get('tsEx')[s1[stack]:s2[stack]]
        tsEy[:,stack] = ts.get('tsEy')[s1[stack]:s2[stack]]
        tsHx[:,stack] = ts.get('tsHx')[s1[stack]:s2[stack]]
        tsHy[:,stack] = ts.get('tsHy')[s1[stack]:s2[stack]]
        tsHz[:,stack] = ts.get('tsHz')[s1[stack]:s2[stack]]
        tsRx[:,stack] = tsR.get('tsRx')[s1[stack]:s2[stack]]
        tsRy[:,stack] = tsR.get('tsRy')[s1[stack]:s2[stack]]
        s1[stack+1] = int(s2[stack]-(WindowLength/2))
        s2[stack+1] = int(s2[stack]+(WindowLength/2))
    if notch == 1:
        print("Applying notch filter...")
        tsEx = notchfil(tsEx,fs)
        tsEy = notchfil(tsEy,fs)
        tsHx = notchfil(tsHx,fs)
        tsHy = notchfil(tsHy,fs)
        tsHz = notchfil(tsHz,fs)
        tsRx = notchfil(tsRx,fs)
        tsRy = notchfil(tsRy,fs)
        print("Done!")
    print("\nApplying detrend...")
    dtsEx = signal.detrend(tsEx,axis= 0)
    dtsEy = signal.detrend(tsEy,axis= 0)
    dtsHx = signal.detrend(tsHx,axis= 0)
    dtsHy = signal.detrend(tsHy,axis= 0)
    dtsHz = signal.detrend(tsHz,axis= 0)
    dtsRx = signal.detrend(tsRx,axis= 0)
    dtsRy = signal.detrend(tsRy,axis= 0)
    print("Done!")
    #
    print("\nComputing FFT...")
    fqs,xfftEx = dofft(dtsEx,fs,WindowLength)
    fqs,xfftEy = dofft(dtsEy,fs,WindowLength)
    fqs,xfftHx = dofft(dtsHx,fs,WindowLength)
    fqs,xfftHy = dofft(dtsHy,fs,WindowLength)
    fqs,xfftHz = dofft(dtsHz,fs,WindowLength)
    fqs,xfftRx = dofft(dtsRx,fs,WindowLength)
    fqs,xfftRy = dofft(dtsRy,fs,WindowLength)
    print("Done!")
    calEx = np.delete(xfftEx, 0, 0)
    calEy = np.delete(xfftEy, 0, 0)
    calHx = np.empty(calEx.shape,dtype=complex)
    calHy = np.empty(calEx.shape,dtype=complex)
    calHz = np.empty(calEx.shape,dtype=complex)
    calRx = np.empty(calEx.shape,dtype=complex)
    calRy = np.empty(calEx.shape,dtype=complex)
    print('\nCalibrating...')
    for stack in range(procinfo.get('nstacks')):
        if ChoppStat == 1:
            f = fqs[:,0]
            # do Calibration
            calHx[:,stack] = calibrateon(f,xfftHx[:,stack],ChoppDataHx,calt_Hx)
            calHy[:,stack] = calibrateon(f,xfftHy[:,stack],ChoppDataHy,calt_Hy)
            calHz[:,stack] = calibrateon(f,xfftHz[:,stack],ChoppDataHz,calt_Hz)
            calRx[:,stack] = calibrateon(f,xfftRx[:,stack],ChoppDataRx,calt_Rx)
            calRy[:,stack] = calibrateon(f,xfftRy[:,stack],ChoppDataRy,calt_Ry)
            f = np.delete(f,0,)        
        elif ChoppStat == 0:
            f = fqs[:,0]
            calHx[:,stack] = calibrateoff(f,xfftHx[:,stack],ChoppDataHx,calt_Hx)
            calHy[:,stack] = calibrateoff(f,xfftHy[:,stack],ChoppDataHy,calt_Hy)
            calHz[:,stack] = calibrateoff(f,xfftHz[:,stack],ChoppDataHz,calt_Hz)
            calRx[:,stack] = calibrateon(f,xfftRx[:,stack],ChoppDataRx,calt_Rx)
            calRy[:,stack] = calibrateon(f,xfftRy[:,stack],ChoppDataRy,calt_Ry)
            f = np.delete(f,0,)
    print("Done!")
    #
    print('\nBand averaging over target frequencies:')
    if config.get('parzen_radius') == 0:
        if fs >= 1000:
            cr = 0.05
        elif (fs < 1000):
            cr = 0.25
        if fs == 512:
            cr = 0.1
        if fs <= 0.125:
            cr = 0.4
    else:
        cr = config.get('parzen_radius')
    for i in tqdm(range(np.shape(ftlist)[0])):
        ft = ftlist[i]
        pf = parzen(f,ft,cr)
        dof[i,0] = (2*2*np.sum(pf!=0))-4
        for stack in (range(procinfo.get('nstacks'))):
            Ex[i,stack] = np.sum((calEx[:,stack]) * pf) / np.sum(pf)
            # Ex[i,stack] = np.sum((calEx * np.conj(calEx)) * pf) / np.sum(pf)
            Ey[i,stack] = np.sum((calEy[:,stack]) * pf) / np.sum(pf)
            # Ey[i,stack] = np.sum((calEy * np.conj(calEy)) * pf) / np.sum(pf)
            Hx[i,stack] = np.sum((calHx[:,stack]) * pf) / np.sum(pf)
            Hy[i,stack] = np.sum((calHy[:,stack]) * pf) / np.sum(pf)
            Hz[i,stack] = np.sum((calHz[:,stack]) * pf) / np.sum(pf)
            Rx[i,stack] = np.sum((calRx[:,stack]) * pf) / np.sum(pf)
            Ry[i,stack] = np.sum((calRy[:,stack]) * pf) / np.sum(pf)
            ExRxc[i,stack] = np.sum((calEx[:,stack] * np.conj(calRx[:,stack])) * pf) / np.sum(pf)
            ExRyc[i,stack] = np.sum((calEx[:,stack] * np.conj(calRy[:,stack])) * pf) / np.sum(pf)
            EyRxc[i,stack] = np.sum((calEy[:,stack] * np.conj(calRx[:,stack])) * pf) / np.sum(pf)
            EyRyc[i,stack] = np.sum((calEy[:,stack] * np.conj(calRy[:,stack])) * pf) / np.sum(pf)
            HxRxc[i,stack] = np.sum((calHx[:,stack] * np.conj(calRx[:,stack])) * pf) / np.sum(pf)
            HxRyc[i,stack] = np.sum((calHx[:,stack] * np.conj(calRy[:,stack])) * pf) / np.sum(pf)
            HyRxc[i,stack] = np.sum((calHy[:,stack] * np.conj(calRx[:,stack])) * pf) / np.sum(pf)
            HyRyc[i,stack] = np.sum((calHy[:,stack] * np.conj(calRy[:,stack])) * pf) / np.sum(pf)
            ExExc[i,stack] = np.sum((calEx[:,stack] * np.conj(calEx[:,stack])) * pf) / np.sum(pf)
            EyEyc[i,stack] = np.sum((calEy[:,stack] * np.conj(calEy[:,stack])) * pf) / np.sum(pf)
            ExEyc[i,stack] = np.sum((calEx[:,stack] * np.conj(calEy[:,stack])) * pf) / np.sum(pf)
            HzHxc[i,stack] = np.sum((calHz[:,stack] * np.conj(calHx[:,stack])) * pf) / np.sum(pf)
            HzHyc[i,stack] = np.sum((calHz[:,stack] * np.conj(calHy[:,stack])) * pf) / np.sum(pf)
            tHxHxc[i,stack] = np.sum((calHx[:,stack] * np.conj(calHx[:,stack])) * pf) / np.sum(pf)
            tHyHyc[i,stack] = np.sum((calHy[:,stack] * np.conj(calHy[:,stack])) * pf) / np.sum(pf)
            tHxHyc[i,stack] = np.sum((calHx[:,stack] * np.conj(calHy[:,stack])) * pf) / np.sum(pf)
            tHyHxc[i,stack] = np.sum((calHy[:,stack] * np.conj(calHx[:,stack])) * pf) / np.sum(pf)
            avgf[i,0] = np.sum(pf!=0)
    Zyy_num = (HxRxc * EyRyc) - (HxRyc * EyRxc)
    Zyx_num = (HyRyc * EyRxc) - (HyRxc * EyRyc)
    Z_deno = (HxRxc * HyRyc) - (HxRyc * HyRxc)
    Zyy_single = -1 * (Zyy_num / Z_deno)
    Zyx_single = -1 * (Zyx_num / Z_deno)
    Zxx_num = (HyRyc * ExRxc) - (HyRxc * ExRyc)
    Zxy_num = (HxRxc * ExRyc) - (HxRyc * ExRxc)
    Zxx_single = -1 * (Zxx_num / Z_deno)
    Zxy_single = -1 * (Zxy_num / Z_deno)
    bandavg = {}
    bandavg['dof'] = dof
    bandavg['avgf'] = avgf
    bandavg['Ex'] = Ex
    bandavg['Ey'] = Ey
    bandavg['Hx'] = Hx
    bandavg['Hy'] = Hy
    bandavg['Hz'] = Hz
    bandavg['ExHxc'] = ExRxc
    bandavg['ExHyc'] = ExRyc
    bandavg['EyHxc'] = EyRxc
    bandavg['EyHyc'] = EyRyc
    bandavg['HxHxc'] = HxRxc
    bandavg['HxHyc'] = HxRyc
    bandavg['HyHxc'] = HyRxc
    bandavg['HyHyc'] = HyRyc
    bandavg['ExExc'] = ExExc
    bandavg['EyEyc'] = EyEyc
    bandavg['ExEyc'] = ExEyc
    bandavg['HzHxc'] = HzHxc
    bandavg['HzHyc'] = HzHyc
    bandavg['tHxHxc'] = tHxHxc
    bandavg['tHyHyc'] = tHyHyc
    bandavg['tHxHxc'] = tHxHxc
    bandavg['tHxHyc'] = tHxHyc
    bandavg['tHyHxc'] = tHyHxc
    bandavg['Zxx_single'] = Zxx_single
    bandavg['Zxy_single'] = Zxy_single
    bandavg['Zyy_single'] = Zyy_single
    bandavg['Zyx_single'] = Zyx_single
    # plt.plot(f,abs(b))
    # plt.xscale('log')
    # plt.gca().invert_xaxis()
    return ftlist,bandavg

def getcalibrationdata(procinfo):
    """
    Function to read calibration data from XML file

    Parameters
    ----------
    procinfo : This dictionary is required to get the path of processing 
                folder

    Returns
    -------
    CalDataHx : Calibaration data for Hx
    CalDataHy : Calibaration data for Hy
    CalDataHz : Calibaration data for Hz

    """
    import xml.etree.ElementTree as ET
    
    xmlfiles = []
    for zz in os.listdir(procinfo['proc_path']):
        if zz.endswith(".xml"):
            xmlfiles.append(zz)
    file = procinfo['proc_path']+'/'+xmlfiles[0]
    mytree = ET.parse(file)
    root = mytree.getroot()
    calsensors = (root[2])
    Hx = calsensors[2]
    Hy = calsensors[3]
    Hz = calsensors[4]
    
    CalDataHx = {}
    freq = []
    for x in Hx.findall(".//caldata[@chopper='on']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hx.findall(".//caldata[@chopper='on']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hx.findall(".//caldata[@chopper='on']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataON = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHx['1'] = dataON
    
    freq = []
    for x in Hx.findall(".//caldata[@chopper='off']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hx.findall(".//caldata[@chopper='off']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hx.findall(".//caldata[@chopper='off']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataOFF = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHx['0'] = dataOFF
    
    CalDataHy = {}
    freq = []
    for x in Hy.findall(".//caldata[@chopper='on']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hy.findall(".//caldata[@chopper='on']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hy.findall(".//caldata[@chopper='on']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataON = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHy['1'] = dataON
    
    freq = []
    for x in Hy.findall(".//caldata[@chopper='off']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hy.findall(".//caldata[@chopper='off']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hy.findall(".//caldata[@chopper='off']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataOFF = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHy['0'] = dataOFF
    
    CalDataHz = {}
    freq = []
    for x in Hz.findall(".//caldata[@chopper='on']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hz.findall(".//caldata[@chopper='on']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hz.findall(".//caldata[@chopper='on']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataON = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHz['1'] = dataON
    
    freq = []
    for x in Hz.findall(".//caldata[@chopper='off']/c1"):
        freq.append(float(x.text))
    
    mag = []
    for x in Hz.findall(".//caldata[@chopper='off']/c2"):
        mag.append(float(x.text))
    
    phase = []
    for x in Hz.findall(".//caldata[@chopper='off']/c3"):
        phase.append(float(x.text))
    
    freq = np.asarray(freq).reshape(-1,1)
    mag = np.asarray(mag).reshape(-1,1)
    phase = np.asarray(phase).reshape(-1,1)
    dataOFF = np.concatenate((freq, mag, phase),axis=1)
    
    CalDataHz['0'] = dataOFF
    
    return CalDataHx, CalDataHy, CalDataHz

def dofft(dts,fs,WindowLength):
    """
    Function to perform FFT.

    Parameters
    ----------
    dts : It is an array of float containing detrended time series data divided
        into a number of windows. Number of sample in a window is the number of
        rows and number of windows is the number of columns.
    fs : It is an integer showing the sampling frequency of measurement.
    WindowLength : It is an integer value representing the window length.

    Returns
    -------
    f : It is an array of float containing frequency list after FFT.
    xfft : It is an array of complex containing data in frequency domain.

    """
    w = np.hanning(WindowLength).reshape(-1,1)
    fft_value = np.fft.fft(dts*w,WindowLength,axis=0)
    xfft = np.asarray(fft_value[0:int(WindowLength/2),:])
    fline = np.asarray([np.linspace(0, int(WindowLength/2), num=int(WindowLength/2),dtype=int)])
    f =  np.asarray(fline * fs/(WindowLength)).T
    return f,xfft

def notchfil(ts,fs):
    """
    Function to perform notch filtering.

    Parameters
    ----------
    ts : It is an array of float containing time series data.
    fs : It is an integer showing the sampling frequency of measurement.

    Returns
    -------
    ts : It is an array of float containing time series data after notch filtering.

    """
    f0 = 50.0  # Frequency to be removed from signal (Hz)
    Qua = 20.0  # Quality factor
    # Design notch filter
    b_notch, a_notch = signal.iirnotch(f0, Qua, fs)
    f1 = 150.0  # Frequency to be removed from signal (Hz)
    # Design notch filter
    b2_notch, a2_notch = signal.iirnotch(f1, Qua, fs)
    f2 = 100.0  # Frequency to be removed from signal (Hz)
    # Design notch filter
    b3_notch, a3_notch = signal.iirnotch(f2, Qua, fs)
    #
    ts = signal.filtfilt(b_notch, a_notch, ts, axis = 0)
    ts = signal.filtfilt(b2_notch, a2_notch, ts, axis = 0)
    ts = signal.filtfilt(b3_notch, a3_notch, ts, axis = 0)
    return ts

def parzen(f,ft,cr):
    """
    Function to create a parzen window

    Parameters
    ----------
    f : It is an array of float containing frequencies.
    ft : It is a float value that is the target frequency.
    cr : It is a float value which is the parzen window radius.

    Returns
    -------
    pf : It is an array of float containing the parzen window values.

    """
    pf = np.empty((np.shape(f)[0],))
    pf[:] = np.nan
    fr = cr * ft
    pf[0] = 0
    for i in range(1, np.shape(f)[0]):
        cond = abs(ft - f[i])
        if cond == 0:
            pf[i] = 1
        elif (cond > 0) and (cond < fr):
            u = (np.pi*abs(cond))/fr
            pf[i] = (np.sin(u)/u) ** 4
        elif (cond > fr) or (cond == fr):
            pf[i] = 0
    # plt.plot(f,pf)
    # plt.xscale('log')
    # plt.gca().invert_xaxis()
    return pf

# Calibration value
def calibrateon(f,xfft,ChoppData,calt):
    """

    Parameters
    ----------
    f : It is an array of float containing frequencies.
    xfft : It is an array of complex containing data in frequency domain.
    ChoppData : It is an array containing calibration information (frequency, 
        magnitude and phase) of the induction coil.
    calt : It is an array of complex containing the calibration values of the induction coil.

    Returns
    -------
    cal : Calibrated spectra.

    """
    minfindx = np.where(f < 0.1)[0]
    cal_all_band = np.interp(f[np.max(minfindx)+1:np.shape(f)[0]],ChoppData[:,0],calt)
    thmag = np.zeros(np.shape(minfindx),)
    thmag[:,] = 0.2 * f[0:np.max(minfindx)+1]
    thph = np.arctan2(4.0,f[0:np.max(minfindx)+1])
    th_band = (thmag * np.cos(thph) + (1j * thmag * np.sin(thph))) * 1000
    cal_all_band = np.concatenate((th_band,cal_all_band))
    xfft = np.delete(xfft,0,)
    cal_all_band = np.delete(cal_all_band,0,)
    cal = xfft/cal_all_band
    return cal

def calibrateoff(f,xfft,ChoppData,calt):
    """

    Parameters
    ----------
    f : It is an array of float containing frequencies.
    xfft : It is an array of complex containing data in frequency domain.
    ChoppData : It is an array containing calibration information (frequency, 
        magnitude and phase) of the induction coil.
    calt : It is an array of complex containing the calibration values of the induction coil.

    Returns
    -------
    cal : Calibrated spectra.

    """
    cal_all_band = np.interp(f[0:np.shape(f)[0]],ChoppData[:,0],calt)
    xfft = np.delete(xfft,0,)
    cal_all_band = np.delete(cal_all_band,0,)
    cal = xfft/cal_all_band
    return cal

def targetfreq(fs):
    """
    It returns target frequencies corresponding to sampling frequency.

    Parameters
    ----------
    fs : It is an integer showing the sampling frequency of measurement.

    Returns
    -------
    ftlist : It is an array of float which is a list of target frequencies.

    """
    if fs == 65536:
        ftlist = ( 1.44850347101214E+0004, 1.11167288155392E+0004, 8.53167852417281E+0003, 
               8.00000000000000E+0003,  5.96636034638832E+0003, 4.44968197286938E+0003, 
               3.31855075962085E+0003,  2.47495870745984E+0003, 1.84581193639211E+0003, 
               1.37659739302252E+0003,  1.02665951233389E+0003, 7.65677575453908E+0002, 
               5.71038540538369E+0002, 4.25877713065949E+0002)
    elif fs == 4096:
        ftlist = ( 7.65677575453908E+0002, 5.71038540538369E+0002, 4.25877713065949E+0002,
                    3.17617487455902E+0002, 2.36877547809548E+0002, 1.76662101025074E+0002,
                    1.31753719283206E+0002, 9.82612707775626E+0001, 7.32827686941220E+0001,
                    5.46539256512696E+0001, 4.07606268475239E+0001, 3.03990734646247E+0001,
                    2.26714783107853E+0001)
        # ftlist = (5.46539256512696E+0001)
    elif fs == 1024:
        ftlist = ( 2.36877547809548E+0002, 1.76662101025074E+0002, 1.31753719283206E+0002,
                  9.82612707775626E+0001,  7.32827686941220E+0001, 5.46539256512696E+0001,
                  4.07606268475239E+0001,  3.03990734646247E+0001, 2.26714783107853E+0001,
                  1.69082761484340E+0001, 1.26101085422250E+0001, 9.40455644624802E+0000,
                  7.01387158203311E+0000, 5.23091066021279E+0000)
    elif fs == 512:
        ftlist = ( 9.82612707775626E+0001, 7.32827686941220E+0001, 5.46539256512696E+0001,
                  4.07606268475239E+0001, 3.03990734646247E+0001, 2.26714783107853E+0001,
                  1.69082761484340E+0001, 1.26101085422250E+0001, 9.40455644624802E+0000,
                  7.01387158203311E+0000, 5.23091066021279E+0000, 3.90118724232419E+0000,
                  2.90948610830488E+0000)
    # elif fs == 256:
    # elif fs == 128:
    elif fs == 32:
        ftlist = ( 7.01387158203311E+0000, 5.23091066021279E+0000, 3.90118724232419E+0000,
                2.90948610830488E+0000, 2.16988031811974E+0000, 1.61828598580477E+0000,
                1.20690966685269E+0000, 9.00107247247826E-0001, 6.71295523434523E-0001,
                5.00648873965966E-0001, 3.73381448636813E-0001, 2.78466033652964E-0001, 
                2.07678587625385E-0001)
    elif fs == 1:
        ftlist =  (2.076785876E-01,  1.548856612E-01,  1.155129584E-01,  8.614899184E-02,
                   6.424949109E-02, 4.791695199E-02,  3.573622528E-02,  2.665189968E-02,
                   1.987685468E-02,  1.482405969E-02,  1.105571024E-02,  8.245293899E-03,
                   6.149299321E-03)
    elif fs==8:
         ftlist = (1.61828598580477E+0000, 1.20690966685269E+0000, 9.00107247247826E-0001,
                   6.71295523434523E-0001, 5.00648873965966E-0001, 3.73381448636813E-0001,
                   2.78466033652964E-0001, 2.07678587625385E-0001, 1.54885661250254E-0001,
                   1.15512958438456E-0001, 8.61489918401509E-0002, 6.42494910995510E-0002,
                   4.79169519964988E-0002)
    elif fs==4:
        ftlist = ( 3.73381448636813E-0001, 2.78466033652964E-0001, 2.07678587625385E-0001,
                   1.54885661250254E-0001, 1.15512958438456E-0001, 8.61489918401509E-0002,
                   6.42494910995510E-0002, 4.79169519964988E-0002, 3.57362252889629E-0002,
                   2.66518996867085E-0002, 1.98768546808371E-0002, 1.48240596973337E-0002,
                   1.10557102438331E-0002)
    elif fs==2:
        ftlist = ( 3.73381448636813E-0001, 2.78466033652964E-0001, 2.07678587625385E-0001,
                   1.54885661250254E-0001, 1.15512958438456E-0001, 8.61489918401509E-0002,
                   6.42494910995510E-0002, 4.79169519964988E-0002, 3.57362252889629E-0002,
                   2.66518996867085E-0002, 1.98768546808371E-0002, 1.48240596973337E-0002,
                   1.10557102438331E-0002)
    elif fs==0.5:
         ftlist = (1.15512958438456E-0001, 8.61489918401509E-0002, 6.42494910995510E-0002,
                   4.79169519964988E-0002, 3.57362252889629E-0002, 2.66518996867085E-0002,
                   1.98768546808371E-0002, 1.48240596973337E-0002, 1.10557102438331E-0002,
                   8.24529389999566E-0003, 6.14929932115646E-0003, 4.58611695347757E-0003,
                   3.42030329189097E-0003, 2.55084524166997E-0003,)
    elif fs==0.125:
         ftlist = (2.66518996867085E-0002, 1.98768546808371E-0002, 1.48240596973337E-0002,
                   1.10557102438331E-0002, 8.24529389999566E-0003, 6.14929932115646E-0003,
                   4.58611695347757E-0003, 3.42030329189097E-0003, 2.55084524166997E-0003,
                   1.90240773745913E-0003, 1.41880626092981E-0003, 1.05813867680238E-0003,
                   7.89154580281698E-0004)
    elif fs==0.03125:
        ftlist = ( 6.14929932115646E-0003, 4.58611695347757E-0003, 3.42030329189097E-0003,
                   2.55084524166997E-0003, 1.90240773745913E-0003, 1.41880626092981E-0003,
                   1.05813867680238E-0003, 7.89154580281698E-0004, 5.88547574370430E-0004,
                   4.38935863710846E-0004)
    elif fs==32768:
        ftlist = ( 8.00000000000000E+0003, 5.96636034638832E+0003, 4.44968197286938E+0003,
         3.31855075962085E+0003, 2.47495870745984E+0003, 1.84581193639211E+0003,
         1.37659739302252E+0003, 1.02665951233389E+0003, 7.65677575453908E+0002,
         5.71038540538369E+0002, 4.25877713065949E+0002, 3.17617487455902E+0002,
         2.36877547809548E+0002)
    elif fs==16384:
        ftlist = ( 3.318550759E+03,  2.474958707E+03,  1.845811936E+03,
                  1.376597393E+03,  1.026659512E+03, 7.656775754E+02,
                  5.710385405E+02,  4.258777130E+02,  3.176174874E+02,
                  2.368775478E+02, 1.766621010E+02,  1.317537192E+02,
                  9.826127077E+01)
    elif fs==128:
        ftlist = ( 3.039907346E+01,  2.267147831E+01,  1.690827614E+01,
                  1.261010854E+01,  9.404556446E+00, 7.013871582E+00,
                  5.230910660E+00,  3.901187242E+00,  2.909486108E+00,
                  2.169880318E+00, 1.618285985E+00,  1.206909666E+00,
                  9.001072472E-01,  6.712955234E-01)
    elif fs==2048:
            ftlist = (1.20000000e+01,  1.50000000e+01,  2.20000000e+01, 
                      3.00000000e+01,  4.00000000e+01, 5.70000000e+01,  
                      7.40000000e+01,  9.80000000e+01,  1.32000000e+02,  
                      1.76000000e+02,  2.36000000e+02,  3.18000000e+02,
                      4.26000000e+02,  5.72000000e+02,  7.66000000e+02)
    return ftlist

def getjackknife(bandavg,mode):
    """
    To get Jackknife values

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values and impedance values for all time windows at all target 
        frequencies.
    mode : It is a string which is either 'Ex' or 'Ey'.

    Returns
    -------
    Array of complex containing Jackknife mean for all target frequencies.

    """
    Z_deno = ((bandavg.get('HxHxc') * bandavg.get('HyHyc')) - 
        ( bandavg.get('HxHyc') * bandavg.get('HyHxc')))
    Zxx_num = ((bandavg.get('HyHyc') * bandavg.get('ExHxc')) - 
    ( bandavg.get('HyHxc') * bandavg.get('ExHyc')))
    Zxy_num = ((bandavg.get('HxHxc') * bandavg.get('ExHyc')) - 
    (bandavg.get('HxHyc') * bandavg.get('ExHxc')))
    Zyx_num = ((bandavg.get('HyHyc') * bandavg.get('EyHxc')) - 
    (bandavg.get('HyHxc') * bandavg.get('EyHyc')))
    Zyy_num = ((bandavg.get('HxHxc') * bandavg.get('EyHyc')) - 
    (bandavg.get('HxHyc') * bandavg.get('EyHxc')))
    Zxx = -1 * (Zxx_num/Z_deno)
    Zxy = -1 * (Zxy_num/Z_deno)
    Zyx = -1 * (Zyx_num/Z_deno)
    Zyy = -1 * (Zyy_num/Z_deno)
    Zxx_jack = np.empty((np.shape(Zxy)[0],1),dtype=complex)
    Zxy_jack = np.empty((np.shape(Zxy)[0],1),dtype=complex)
    Zyx_jack = np.empty((np.shape(Zxy)[0],1),dtype=complex)
    Zyy_jack = np.empty((np.shape(Zxy)[0],1),dtype=complex)
    if (mode == 'Ex'):
        for i in range(np.shape(Zxy)[0]):
            Zxx_jack[i,0] = jackknife(Zxx[i,:])
            Zxy_jack[i,0] = jackknife(Zxy[i,:])
        return Zxx_jack,Zxy_jack
    elif (mode == 'Ey'):
        for i in range(np.shape(Zxy)[0]):
            Zyx_jack[i,0] = jackknife(Zyx[i,:])
            Zyy_jack[i,0] = jackknife(Zyy[i,:])
        return Zyx_jack,Zyy_jack

def jackknife(Z):
    """

    Parameters
    ----------
    Z : Array of complex which are the impedance values for a target frequency.

    Returns
    -------
    mean_jackknife : Complex value which is the Jackknife mean for a target frequency.

    """
    nstacks = np.shape(Z)[0]
    for k in range(nstacks-1):
        Zminusi = np.empty((np.shape(Z)[0],),dtype=complex)
        Zidiff = np.empty((np.shape(Z)[0],),dtype=complex)
        for j in range(np.shape(Z)[0]):
            Zminusi[j] = (np.sum(Z)-Z[j])/(np.shape(Z)[0]-1)
        for j in range(np.shape(Z)[0]):
            Zidiff[j] = abs(Zminusi[j] - np.mean(Zminusi))
        mean_jackknife = (np.shape(Z)[0] * np.mean(Z)) - (((np.shape(Z)[0]-1)/np.shape(Z)[0])*np.sum(Zminusi))
        Zvar = ((np.shape(Z)[0]-1)/(np.shape(Z)[0]) * np.sum(Zidiff ** 2))
        ind = (np.where(Zidiff == np.max(Zidiff))[0])[0]
        Z = np.delete(Z,ind)
        del Zminusi, Zidiff
    return mean_jackknife

def huberEx(bandavg,Z_jackk,stacki):
    """
    Computing huber estimates for Ex.

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.
    Z_jackk : Array of complex containing Jackknife mean for all target frequencies.
    stacki : It is an integer value to select all windows for a target frequency.

    Returns
    -------
    Zxx_robust_huber : It is a complex value which is the huber mean for a target frequency.
    Zxy_robust_huber : It is a complex value which is the huber mean for a target frequency.
    bandavgEx_huber : It is a Python dictionary containing averaged value of auto and cross
        spectra.

    """
    Ex = bandavg.get('Ex')
    Hx = bandavg.get('Hx')
    Hy = bandavg.get('Hy')
    nstacks = Hx.shape[1]
    Zxx_jackk = Z_jackk.get('Zxx')[stacki,0]
    Zxy_jackk = Z_jackk.get('Zxy')[stacki,0]
    rxl = abs(abs(Ex) - (abs(Zxx_jackk) * abs(Hx)) - (abs(Zxy_jackk) * abs(Hy)))
    dmx = 1.483 * np.median(abs(rxl - np.median(rxl,axis=1).reshape(-1, 1)),axis=1).reshape(-1, 1)
    kmx = 1.5 * dmx
    huber_matrix = huberwt(rxl,kmx)
    ExExc_hup = bandavg.get('ExExc') * huber_matrix
    HxHxc_hup = bandavg.get('HxHxc') * huber_matrix
    ExHyc_hup = bandavg.get('ExHyc') * huber_matrix
    ExHxc_hup = bandavg.get('ExHxc') * huber_matrix
    HxHyc_hup = bandavg.get('HxHyc') * huber_matrix
    HyHyc_hup = bandavg.get('HyHyc') * huber_matrix
    HyHxc_hup = bandavg.get('HyHxc') * huber_matrix
    ExExc_hup_avg = (np.sum(ExExc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HxHxc_hup_avg = (np.sum(HxHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    ExHyc_hup_avg = (np.sum(ExHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    ExHxc_hup_avg = (np.sum(ExHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HxHyc_hup_avg = (np.sum(HxHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HyHyc_hup_avg = (np.sum(HyHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HyHxc_hup_avg = (np.sum(HyHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    Zxx_new = (HyHyc_hup_avg * ExHxc_hup_avg) - (HyHxc_hup_avg * ExHyc_hup_avg)
    Zxy_new = (HxHxc_hup_avg * ExHyc_hup_avg) - (HxHyc_hup_avg * ExHxc_hup_avg)
    Z_deno_new = (HxHxc_hup_avg * HyHyc_hup_avg) - (HxHyc_hup_avg * HyHxc_hup_avg)
    Zxx_robust_huber = -1 * (Zxx_new / Z_deno_new)
    Zxy_robust_huber = -1 * (Zxy_new / Z_deno_new)
    for i in range(4):
        Lc = np.sum((huber_matrix==1)*1,axis=1)
        for kk in range(np.size(Lc)):
            if Lc[kk] == 0:
                Lc[kk] = 1
        dhx = (np.sqrt((nstacks/(Lc ** 2))*(np.sum(huber_matrix * (rxl ** 2),axis=1)))).reshape(-1,1)
        khx = 1.5 * dhx
        Ex_hup = Ex * huber_matrix
        Hx_hup = Hx * huber_matrix
        Hy_hup = Hy * huber_matrix
        rxl = abs(abs(Ex_hup) - (abs(Zxx_robust_huber) * abs(Hx_hup)) - (abs(Zxy_robust_huber) * abs(Hy_hup)))
        huber_matrix = huberwt(rxl,khx)
        ExExc_hup = ExExc_hup * huber_matrix
        HxHxc_hup = HxHxc_hup * huber_matrix
        ExHyc_hup = ExHyc_hup * huber_matrix
        ExHxc_hup = ExHxc_hup * huber_matrix
        HxHyc_hup = HxHyc_hup * huber_matrix
        HyHyc_hup = HyHyc_hup * huber_matrix
        HyHxc_hup = HyHxc_hup * huber_matrix
        ExExc_hup_avg = np.sum(ExExc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HxHxc_hup_avg = np.sum(HxHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        ExHyc_hup_avg = np.sum(ExHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        ExHxc_hup_avg = np.sum(ExHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HxHyc_hup_avg = np.sum(HxHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HyHyc_hup_avg = np.sum(HyHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HyHxc_hup_avg = np.sum(HyHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        Zxx_new = (HyHyc_hup_avg * ExHxc_hup_avg) - (HyHxc_hup_avg * ExHyc_hup_avg)
        Zxy_new = (HxHxc_hup_avg * ExHyc_hup_avg) - (HxHyc_hup_avg * ExHxc_hup_avg)
        Z_deno_new = (HxHxc_hup_avg * HyHyc_hup_avg) - (HxHyc_hup_avg * HyHxc_hup_avg)
        Zxx_robust_huber = -1 * (Zxx_new / Z_deno_new).reshape(-1,1)
        Zxy_robust_huber = -1 * (Zxy_new / Z_deno_new).reshape(-1,1)
    bandavgEx_huber = {}
    bandavgEx_huber['rxl'] = rxl
    bandavgEx_huber['khx'] = khx
    bandavgEx_huber['ExExc'] = ExExc_hup_avg
    bandavgEx_huber['HxHxc'] = HxHxc_hup_avg
    bandavgEx_huber['ExHyc'] = ExHyc_hup_avg
    bandavgEx_huber['ExHxc'] = ExHxc_hup_avg
    bandavgEx_huber['HxHyc'] = HxHyc_hup_avg
    bandavgEx_huber['HyHyc'] = HyHyc_hup_avg
    bandavgEx_huber['HyHxc'] = HyHxc_hup_avg
    bandavgEx_huber['Ex'] = Ex_hup
    bandavgEx_huber['huber_matrix'] = huber_matrix
    return Zxx_robust_huber, Zxy_robust_huber,bandavgEx_huber

def huberEy(bandavg,Z_jackk,stacki):
    """
    Computing huber estimates for Ey.

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.
    Z_jackk : Array of complex containing Jackknife mean for all target frequencies.
    stacki : It is an integer value to select all windows for a target frequency.

    Returns
    -------
    Zyy_robust_huber : It is a complex value which is the huber mean for a target frequency.
    Zyx_robust_huber : It is a complex value which is the huber mean for a target frequency.
    bandavgEy_huber : It is a Python dictionary containing averaged value of auto and cross
        spectra.

    """
    Ey = bandavg.get('Ey')
    Hx = bandavg.get('Hx')
    Hy = bandavg.get('Hy')
    nstacks = Hx.shape[1]
    Zyy_jackk = Z_jackk.get('Zyy')[stacki,0]
    Zyx_jackk = Z_jackk.get('Zyx')[stacki,0]
    ryl = abs(abs(Ey) - (abs(Zyy_jackk) * abs(Hy)) - (abs(Zyx_jackk) * abs(Hx)))
    dmy = 1.483 * np.median(abs(ryl - np.median(ryl,axis=1).reshape(-1, 1)),axis=1).reshape(-1, 1)
    kmy = 1.5 * dmy
    huber_matrix = huberwt(ryl,kmy)
    EyEyc_hup = bandavg.get('EyEyc') * huber_matrix
    HxHxc_hup = bandavg.get('HxHxc') * huber_matrix
    EyHyc_hup = bandavg.get('EyHyc') * huber_matrix
    EyHxc_hup = bandavg.get('EyHxc') * huber_matrix
    HxHyc_hup = bandavg.get('HxHyc') * huber_matrix
    HyHyc_hup = bandavg.get('HyHyc') * huber_matrix
    HyHxc_hup = bandavg.get('HyHxc') * huber_matrix
    EyEyc_hup_avg = (np.sum(EyEyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HxHxc_hup_avg = (np.sum(HxHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    EyHyc_hup_avg = (np.sum(EyHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    EyHxc_hup_avg = (np.sum(EyHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HxHyc_hup_avg = (np.sum(HxHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HyHyc_hup_avg = (np.sum(HyHyc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    HyHxc_hup_avg = (np.sum(HyHxc_hup,axis=1)/np.sum(huber_matrix,axis=1)).reshape(-1,1)
    Zyy_new = (HxHxc_hup_avg * EyHyc_hup_avg) - (HxHyc_hup_avg * EyHxc_hup_avg)
    Zyx_new = (HyHyc_hup_avg * EyHxc_hup_avg) - (HyHxc_hup_avg * EyHyc_hup_avg)
    Z_deno_new = (HxHxc_hup_avg * HyHyc_hup_avg) - (HxHyc_hup_avg * HyHxc_hup_avg)
    Zyy_robust_huber = -1 * (Zyy_new / Z_deno_new)
    Zyx_robust_huber = -1 * (Zyx_new / Z_deno_new)
    for i in range(4):
        Lc = np.sum((huber_matrix==1)*1,axis=1)
        for kk in range(np.size(Lc)):
            if Lc[kk] == 0:
                Lc[kk] = 1
        dhy = (np.sqrt((nstacks/(Lc ** 2))*(np.sum(huber_matrix * (ryl ** 2),axis=1)))).reshape(-1,1)
        khy = 1.5 * dhy
        Ey_hup = Ey * huber_matrix
        Hx_hup = Hx * huber_matrix
        Hy_hup = Hy * huber_matrix
        ryl = abs(abs(Ey_hup) - (abs(Zyy_robust_huber) * abs(Hy_hup)) - (abs(Zyx_robust_huber) * abs(Hx_hup)))
        huber_matrix = huberwt(ryl,khy)
        EyEyc_hup = EyEyc_hup * huber_matrix
        HxHxc_hup = HxHxc_hup * huber_matrix
        EyHyc_hup = EyHyc_hup * huber_matrix
        EyHxc_hup = EyHxc_hup * huber_matrix
        HxHyc_hup = HxHyc_hup * huber_matrix
        HyHyc_hup = HyHyc_hup * huber_matrix
        HyHxc_hup = HyHxc_hup * huber_matrix
        EyEyc_hup_avg = np.sum(EyEyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HxHxc_hup_avg = np.sum(HxHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        EyHyc_hup_avg = np.sum(EyHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        EyHxc_hup_avg = np.sum(EyHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HxHyc_hup_avg = np.sum(HxHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HyHyc_hup_avg = np.sum(HyHyc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        HyHxc_hup_avg = np.sum(HyHxc_hup,axis=1) / np.sum(huber_matrix,axis=1)
        Zyy_new = (HxHxc_hup_avg * EyHyc_hup_avg) - (HxHyc_hup_avg * EyHxc_hup_avg)
        Zyx_new = (HyHyc_hup_avg * EyHxc_hup_avg) - (HyHxc_hup_avg * EyHyc_hup_avg)
        Z_deno_new = (HxHxc_hup_avg * HyHyc_hup_avg) - (HxHyc_hup_avg * HyHxc_hup_avg)
        Zyy_robust_huber = -1 * (Zyy_new / Z_deno_new).reshape(-1,1)
        Zyx_robust_huber = -1 * (Zyx_new / Z_deno_new).reshape(-1,1)
    bandavgEy_huber = {}
    bandavgEy_huber['ryl'] = ryl
    bandavgEy_huber['khy'] = khy
    bandavgEy_huber['EyEyc'] = EyEyc_hup_avg
    bandavgEy_huber['HxHxc'] = HxHxc_hup_avg
    bandavgEy_huber['EyHyc'] = EyHyc_hup_avg
    bandavgEy_huber['EyHxc'] = EyHxc_hup_avg
    bandavgEy_huber['HxHyc'] = HxHyc_hup_avg
    bandavgEy_huber['HyHyc'] = HyHyc_hup_avg
    bandavgEy_huber['HyHxc'] = HyHxc_hup_avg
    bandavgEy_huber['huber_matrix'] = huber_matrix
    return Zyy_robust_huber, Zyx_robust_huber,bandavgEy_huber

def huberwt(rl,km):
    """

    Parameters
    ----------
    rl : It is an array of float which are the residuals. Number of rows are number
        of target frequencies and number of columns are number of time windows.
    km : It is an array of float. Number of rows are number of target frequencies.

    Returns
    -------
    huber_matrix : It is an array of float containing huber weights. Number of 
    rows are number of target frequencies and number of columns are number of time windows.

    """
    huber_matrix1 = (rl <= km) * 1
    huber_matrix2 = (rl > km) * 1
    huber_matrix2 = huber_matrix2 * (km/rl)
    huber_matrix = huber_matrix1 + huber_matrix2
    return huber_matrix


def perform_robust(ftlist,bandavg):
    """

    Parameters
    ----------
    ftlist : It is an array of float which is a list of target frequencies.
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    Z_huber : It is a Python dictionary containing huber estimates for all target
        frequencies.

    """
    Zxx_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zxy_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zyx_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zyy_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    for stacki in range(np.shape(ftlist)[0]):
        print('\nComputing Jackknife estimate....ft = ' + str(ftlist[stacki,0]))
        bandavg_singleEx = makeband(bandavg,stacki,'selectedEx')
        bandavg_singleEy = makeband(bandavg,stacki,'selectedEy')
        # bandavgT_single = tipper.makeband(bandavg,stacki)
        Zxx_jackk[stacki,0],Zxy_jackk[stacki,0] = (getjackknife(bandavg_singleEx,'Ex'))
        Zyx_jackk[stacki,0],Zyy_jackk[stacki,0] = (getjackknife(bandavg_singleEy,'Ey'))
        # Tx_jackk[stacki,0],Ty_jackk[stacki,0] = (tipper.getjackknife(bandavgT_single))
    print('Finished.')
    Z_jackk = {'Zxx': Zxx_jackk}
    Z_jackk['Zxy'] = Zxy_jackk
    Z_jackk['Zyx'] = Zyx_jackk
    Z_jackk['Zyy'] = Zyy_jackk
    del Zxx_jackk,Zxy_jackk,Zyx_jackk,Zyy_jackk
    #
    Z_huber = {}
    Zxx_huber = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zxy_huber = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zyx_huber = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    Zyy_huber = np.empty((np.shape(ftlist)[0],1),dtype=complex)
    for stacki in range(np.shape(ftlist)[0]):
        print('\nComputing Huber estimate....ft = ' + str(ftlist[stacki,0]))
        bandavg_singleEx = makeband(bandavg,stacki,'selectedEx')
        bandavg_singleEy = makeband(bandavg,stacki,'selectedEy')
        Zxx_huber[stacki,0],Zxy_huber[stacki,0],bandavgEx_huber = huberEx(bandavg_singleEx,Z_jackk,stacki)
        Zyy_huber[stacki,0],Zyx_huber[stacki,0],bandavgEy_huber = huberEy(bandavg_singleEy,Z_jackk,stacki)
    Z_huber['Zxx'] = Zxx_huber
    Z_huber['Zxy'] = Zxy_huber
    Z_huber['Zyy'] = Zyy_huber
    Z_huber['Zyx'] = Zyx_huber
    return Z_huber

def makeband(bandavg,i,coh_mode):
    """

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.
    i : It is an integer to select all windows for a target frequency.
    coh_mode : It is a string value.

    Returns
    -------
    bandavg_single : It is a Python dictionary containing selected auto and cross spectra
        values for a target frequency.

    """
    Ex = bandavg.get('Ex')[i,:].reshape(1,-1)
    Ey = bandavg.get('Ey')[i,:].reshape(1,-1)
    Hx = bandavg.get('Hx')[i,:].reshape(1,-1)
    Hy = bandavg.get('Hy')[i,:].reshape(1,-1)
    Hz = bandavg.get('Hz')[i,:].reshape(1,-1)
    ExHxc = bandavg.get('ExHxc')[i,:].reshape(1,-1)
    ExHyc = bandavg.get('ExHyc')[i,:].reshape(1,-1)
    EyHxc = bandavg.get('EyHxc')[i,:].reshape(1,-1)
    EyHyc = bandavg.get('EyHyc')[i,:].reshape(1,-1)
    HxHxc = bandavg.get('HxHxc')[i,:].reshape(1,-1)
    HxHyc = bandavg.get('HxHyc')[i,:].reshape(1,-1)
    HyHxc = bandavg.get('HyHxc')[i,:].reshape(1,-1)
    HyHyc = bandavg.get('HyHyc')[i,:].reshape(1,-1)
    ExExc = bandavg.get('ExExc')[i,:].reshape(1,-1)
    EyEyc = bandavg.get('EyEyc')[i,:].reshape(1,-1)
    ExEyc = bandavg.get('ExEyc')[i,:].reshape(1,-1)
    coh_selected = bandavg.get(coh_mode)[i,:].reshape(1,-1)
    ind_coh = np.where(coh_selected==0)[1].reshape(1,-1)
    Ex = np.delete(Ex,ind_coh).reshape(1,-1)
    Ey = np.delete(Ey,ind_coh).reshape(1,-1)
    Hx = np.delete(Hx,ind_coh).reshape(1,-1)
    Hy = np.delete(Hy,ind_coh).reshape(1,-1)
    Hz = np.delete(Hz,ind_coh).reshape(1,-1)
    ExHxc = np.delete(ExHxc,ind_coh).reshape(1,-1)
    ExHyc = np.delete(ExHyc,ind_coh).reshape(1,-1)
    EyHxc = np.delete(EyHxc,ind_coh).reshape(1,-1)
    EyHyc = np.delete(EyHyc,ind_coh).reshape(1,-1)
    HxHxc = np.delete(HxHxc,ind_coh).reshape(1,-1)
    HxHyc = np.delete(HxHyc,ind_coh).reshape(1,-1)
    HyHxc = np.delete(HyHxc,ind_coh).reshape(1,-1)
    HyHyc = np.delete(HyHyc,ind_coh).reshape(1,-1)
    ExExc = np.delete(ExExc,ind_coh).reshape(1,-1)
    EyEyc = np.delete(EyEyc,ind_coh).reshape(1,-1)
    ExEyc = np.delete(ExEyc,ind_coh).reshape(1,-1)
    bandavg_single = {}
    bandavg_single['Ex'] = Ex
    bandavg_single['Ey'] = Ey
    bandavg_single['Hx'] = Hx
    bandavg_single['Hy'] = Hy
    bandavg_single['Hz'] = Hz
    bandavg_single['ExHxc'] = ExHxc
    bandavg_single['ExHyc'] = ExHyc
    bandavg_single['EyHxc'] = EyHxc
    bandavg_single['EyHyc'] = EyHyc
    bandavg_single['HxHxc'] = HxHxc
    bandavg_single['HxHyc'] = HxHyc
    bandavg_single['HyHxc'] = HyHxc
    bandavg_single['HyHyc'] = HyHyc
    bandavg_single['ExExc'] = ExExc
    bandavg_single['EyEyc'] = EyEyc
    bandavg_single['ExEyc'] = ExEyc
    return bandavg_single
    
    
def measid(siteindex):
    """

    Parameters
    ----------
    siteindex : It is an integer representing the serial number of selected site.

    Returns
    -------
    measid : It is a Python dictionary containing measurement ID.

    """
    measid = {}
    measid['Hx'] = 10 * siteindex + 3 + (1/1000)
    measid['Hy'] = 10 * siteindex + 4 + (1/1000)
    measid['Hz'] = 10 * siteindex + 5 + (1/1000)
    measid['Ex'] = 10 * siteindex + 1 + (1/1000)
    measid['Ey'] = 10 * siteindex + 2 + (1/1000)
    return measid

def cleanSpec(bandavg):
    """

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    spmat : It is an array of float containing zeros and ones.

    """
    HxHxc = bandavg.get('HxHxc')
    HyHyc = bandavg.get('HyHyc')
    HxHxc_real = HxHxc.real
    HyHyc_real = HyHyc.real
    HxHxc_imag = HxHxc.imag
    HyHyc_imag = HyHyc.imag
    spmatHx = np.ones(HxHxc.shape)
    spmatHy = np.ones(HyHyc.shape)
    spmatHxi = np.ones(HxHxc.shape)
    spmatHyi = np.ones(HyHyc.shape)
    for i in range(np.shape(HxHxc)[0]):
        for j in range(np.shape(HxHxc)[1]):
            if((HxHxc_real[i,j] < 0) or (HxHxc_real[i,j] >= 5000)):
                spmatHx[i,j] = 0
            if((HyHyc_real[i,j] < 0) or (HyHyc_real[i,j] >= 5000)):
                spmatHy[i,j] = 0
    for i in range(np.shape(HxHxc)[0]):
        for j in range(np.shape(HxHxc)[1]):
            if(HxHxc_imag[i,j] < -100):
                spmatHxi[i,j] = 0
            if(HyHyc_imag[i,j] < -100):
                spmatHyi[i,j] = 0
    spmat = spmatHx * spmatHy * spmatHxi * spmatHyi
    return spmat