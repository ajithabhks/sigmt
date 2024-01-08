# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:27:42 2021
@author: ajithabh

##########
Program: SigMT
Program for the processing of MT time series.
Authors: Ajithabh K.S. and Prasanta K. Patro
CSIR - National Geophysical Research Institute, Hyderabad, India.
##########

This code controls remote reference method in MT data processing.

Give the path of the folder where your sites are kept in the
'project_path' variable.

Then inputs such as site name and measurements will be asked in the console.
First, the details of the local site will be asked. Then, the details of the 
remote site need to be entered.

In case decimation is required, make 'dflag = 1' and give decimation sequence 
in 'decimate' variable.

Make 'ctflag =1' to enable coherency threshold.

Make 'pdflag =1' to enable polarization direction based data selection.

Read user manual for more details about the data selection
"""
# Importing necessary modules
import mtprocRR, var, data_sel, tipper, mahaDist, mtproc, config, plotting
from matplotlib import pyplot as plt
from scipy import signal
import numpy as np
import os
import time
import math
config = config.configuration()
# define project path where sites are kept
project_path = 'D:/NGRI/FIELD RAW DATA/demo/'
# #
# #========= Selection of site and setting a path =========
sites, selectedsite, measid, all_meas, select_meas, proc_path = mtproc.makeprocpath(project_path)
# #========= Site is selected and path is created =========
# #========= Selection of remote site and setting a path ==
sites, selectedsiteR, measidR, all_measR, select_measR, proc_pathR = mtproc.makeprocpath(project_path)
# #========= Remote site is selected and path is created ==
#-------- Time series reading starts ------------------
#
#
procinfo = {}
procinfoR = {}
ts, procinfo['fs'], procinfo['sensor_no'], timeline, procinfo['ChoppStat'], loc = mtprocRR.ts(proc_path)
tsR, procinfoR['fs'], procinfoR['sensor_no'], timelineR, procinfoR['ChoppStat'], locR = mtprocRR.ts(proc_pathR)
#==============
#==============
# Prepare ts data for RR 
timeline = np.asarray(timeline)
timelineR = np.asarray(timelineR)
del tsR['tsEx']
del tsR['tsEy']
del tsR['tsHz']
tsRx = tsR.get('tsHx')
tsRy = tsR.get('tsHy')
Rstart = timeline[0]
Rend = timeline[-1]
Rstart_ind = np.where(np.equal(timelineR,Rstart))[0][0]
Rend_ind = np.where(np.equal(timelineR,Rend))[0][0]
tsRx = tsRx[Rstart_ind:Rend_ind+1]
tsRy = tsRy[Rstart_ind:Rend_ind+1]
del tsR['tsHx']
del tsR['tsHy']
tsR['tsRx'] = tsRx
tsR['tsRy'] = tsRy 
del timelineR, Rstart, Rend
del Rstart_ind, Rend_ind, tsRx, tsRy
#========= Decimation section ================= 
# Keep dflag = 0 if decimation is not required
dflag = 0
if dflag == 1:
    decimate = [8,8,4]
    for d in decimate:
        ts['tsEx'] = signal.decimate(ts.get('tsEx'), d, n=None, ftype='iir')
        ts['tsEy'] = signal.decimate(ts.get('tsEy'), d, n=None, ftype='iir')
        ts['tsHx'] = signal.decimate(ts.get('tsHx'), d, n=None, ftype='iir')
        ts['tsHy'] = signal.decimate(ts.get('tsHy'), d, n=None, ftype='iir')
        ts['tsHz'] = signal.decimate(ts.get('tsHz'), d, n=None, ftype='iir')
        tsR['tsRx'] = signal.decimate(tsR.get('tsRx'), d, n=None, ftype='iir')
        tsR['tsRy'] = signal.decimate(tsR.get('tsRy'), d, n=None, ftype='iir')
        procinfo['fs'] = procinfo.get('fs')/d
#========= Decimation section end =============
#
# Some calculations and printing some information
# No need to edit
procinfo['nofs'] = len(ts['tsEx'])
procinfo['notch'] = 0 # Notch flag 1 - On, 0 - Off
procinfo['notch_freqs'] = [50,100] # Provide frequencies to be filtered seperted by a comma, Eg:[50, 100]
print('--------------------')
#print('\nSaved time series after trend & bias removal.')
print('MT site: ' + selectedsite)
print('Measurement directory: ' + all_meas[select_meas])
print('fs= '+ str(procinfo.get('fs')) +' Hz')
print('Sensor numbers:')
print(procinfo.get('sensor_no'))
print('Length of time series = ' + str(procinfo.get('nofs')))
print('--------------------')
procinfo['meas'] = all_meas[select_meas]
procinfo['proc_path'] = proc_path
procinfoR['proc_path'] = proc_pathR
procinfo['selectedsite'] = selectedsite
del all_meas, select_meas, selectedsite
del all_measR, select_measR, selectedsiteR, proc_pathR
del proc_path, project_path, sites
print('Unused variables deleted.')
print('\n\n--------------------')
# Find out window length
if config.get('FFT_Length') == 0:
    procinfo['WindowLength'] = mtproc.FFTLength(procinfo.get('nofs'))
else:
    procinfo['WindowLength'] = config.get('FFT_Length')
procinfo['overlap'] = 50 #Input in percentage 
print('\nWindow Length selected: '+ str(procinfo.get('WindowLength')))
procinfo['nstacks'] = math.floor(procinfo.get('nofs')/procinfo.get('WindowLength'))
procinfo['nstacks'] = (procinfo.get('nstacks') * 2) - 1
print('Time series overlap: ' + str(procinfo.get('overlap'))+'%')
print('No. of windows: '+ str(procinfo.get('nstacks')))
print('--------------------')
#==================== Start band averaging ====================
# No need to edit
# Band average value after calibration and averaging using parzen window
ftlist,bandavg = mtprocRR.bandavg(ts,procinfo,tsR,procinfoR,config)
#
#==================== Band averaging finished ====================
#########
spmat = mtprocRR.cleanSpec(bandavg)
#
#====Data selection tools section. Coherency threshold & Polarization direction
# Calculation of coherency values for all time windows
AllcohEx = data_sel.cohEx(bandavg)
AllcohEy = data_sel.cohEy(bandavg)
# Calculation of polarization directions for all time windows
alpha_degH,alpha_degE = data_sel.pdvalues(bandavg)
#
#====== Coherency threshold ======
ctflag = 0 # Give '1' to perform coherency threshold based selection
minpercent = 20 # Minimum percentage of windows required
CohThre = 0.9 # Coherency Threshold value Range: (0,1)
[cohMatrixEx, cohMatrixEy] = data_sel.performct(ctflag,CohThre,minpercent,ftlist,bandavg,AllcohEx,AllcohEy)
#
#====== Polarization direction ======
pdflag = 0 # Give '1' to perform polarization direction based selection
pdlim = [-10,10] # Ploarization direction limit
alpha = alpha_degE # Use either alpha_degE (electric field) or alpha_degH (magnetic field)
pdmat = data_sel.performpd(pdflag,pdlim,alpha,bandavg)

# This can be used to mask time windows based on the polarization directions
mwflag = 0 # Give '1' to perform mask windows
timewindow_limits = [0,40] #Time window limits
mwmat = data_sel.performmw(mwflag,timewindow_limits,bandavg)

#====End of data selection tools section
#
bandavg['cohMatrixEx'] = cohMatrixEx
bandavg['cohMatrixEy'] = cohMatrixEy
bandavg['pre_sel_matEx'] = cohMatrixEx * pdmat
bandavg['pre_sel_matEy'] = cohMatrixEy * pdmat
bandavg['mdmatrixEx'],bandavg['Zxx_mcd'],bandavg['Zxy_mcd'],bandavg['mahal_robustEx'] = mahaDist.mcd(bandavg,'Ex',config)
bandavg['mdmatrixEy'],bandavg['Zyx_mcd'],bandavg['Zyy_mcd'],bandavg['mahal_robustEy'] = mahaDist.mcd(bandavg,'Ey',config)
bandavg['selectedEx'] = bandavg.get('mdmatrixEx') * spmat
bandavg['selectedEy'] = bandavg.get('mdmatrixEy') * spmat

# bandavg['coh_selected'] = selmatEx * selmatEy
bandavg['avgt'] = np.sum((bandavg.get('selectedEx'))!=0,axis=1)
del cohMatrixEx,cohMatrixEy
#
#=========== Tipper estimation ===================================
[TxAll, TyAll] = tipper.tippall(bandavg)
mahaWtTx, Tx_mcd_mean = tipper.mcd(TxAll,config)
mahaWtTy, Ty_mcd_mean = tipper.mcd(TyAll,config)
bandavg['tipp_selected'] = mahaWtTx * mahaWtTy
[Tx, Ty] = tipper.tipper(bandavg)
[TxVar, TyVar] = tipper.tipperVar(bandavg)
#=========== Tipper estimation DONE===============================
#==================== Robust estimation begins ===================
#
Z_huber = mtprocRR.perform_robust(ftlist,bandavg)
#
#==================== Calculation of variance ===================
Zvar = {}
Zvar['xx'],Zvar['xy'],cohEx = var.ZExvar(Z_huber,bandavg)
Zvar['yx'],Zvar['yy'],cohEy = var.ZEyvar(Z_huber,bandavg)
# ----------------------------------------------------------------------
#
#
### Plotting figures ###
plotting.plotfigs(procinfo, ftlist, Z_huber, Zvar, Tx, Ty, cohEx, cohEy)
print('Finished.')