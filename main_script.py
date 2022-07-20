# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:49:35 2020

@author: AJITHABH
"""
import mtproc, var, coherency, tipper, mahaDist
from scipy import signal
from matplotlib import pyplot as plt
import numpy as np
import os
import time
import math
#
#
# provide project path where sites are kept
#
project_path = 'D:/NGRI/FIELD RAW DATA/demo/'
#
#
os.chdir(project_path)
sites = [d for d in os.listdir('.') if os.path.isdir(d)]
# all site names are stored in variable 'sites'
print('Sites in the project are:  ')
print(sites)
selectedsite = input("Enter the site: ")
siteindex = sites.index(selectedsite)
measid = mtproc.measid(siteindex)
del siteindex
os.chdir(project_path+selectedsite)
all_meas = [d for d in os.listdir('.') if os.path.isdir(d)]
for i in range(len(all_meas)):
    '\n'
    print ([list((i, all_meas[i]))]) 
select_meas = input("Select measurement: ")
select_meas = int(select_meas)
proc_path = project_path+selectedsite+'/'+all_meas[select_meas]
#
#
# read all channel from the selected measurement
timer_start = time.time()
#
#
#-------- Time series reading starts ------------------
#
#
procinfo = {}
[ts,procinfo['fs'],procinfo['sensor_no'],timeline,procinfo['ChoppStat'],loc] = mtproc.ts(proc_path)
# trend removed 
# V - mv/km
# H - mV
dflag = 0
if dflag == 1:
    decimate = [8,4]
    for d in decimate:
        ts['tsEx'] = signal.decimate(ts.get('tsEx'), d, n=None, ftype='iir')
        ts['tsEy'] = signal.decimate(ts.get('tsEy'), d, n=None, ftype='iir')
        ts['tsHx'] = signal.decimate(ts.get('tsHx'), d, n=None, ftype='iir')
        ts['tsHy'] = signal.decimate(ts.get('tsHy'), d, n=None, ftype='iir')
        ts['tsHz'] = signal.decimate(ts.get('tsHz'), d, n=None, ftype='iir')
        procinfo['fs'] = procinfo.get('fs')/d
procinfo['fs'] = procinfo.get('fs')
procinfo['nofs'] = len(ts['tsEx'])
procinfo['notch'] = 0 # Notch flag 1 - On, 0 - Off
print('--------------------')
print('MT site: ' + selectedsite)
print('Measurement directory: ' + all_meas[select_meas])
print('fs= '+ str(procinfo.get('fs')) +' Hz')
print('Sensor numbers:')
print(procinfo.get('sensor_no'))
print('Length of time series = ' + str(procinfo.get('nofs')))
print('--------------------')
procinfo['meas'] = all_meas[select_meas]
procinfo['proc_path'] = proc_path
procinfo['selectedsite'] = selectedsite
del all_meas, i, select_meas, selectedsite
del proc_path, project_path, sites
print('Unused variables deleted.')
print('\n\n--------------------')
# Find out window length
procinfo['WindowLength'] = mtproc.FFTLength(procinfo.get('nofs'))
#procinfo['WindowLength'] = 512
procinfo['overlap'] = 50 #Input in percentage %
print('\nWindow Length selected: '+ str(procinfo.get('WindowLength')))
procinfo['nstacks'] = math.floor(procinfo.get('nofs')/procinfo.get('WindowLength'))
procinfo['nstacks'] = (procinfo.get('nstacks') * 2) - 1
# procinfo['nstacks'] =300
# if procinfo.get('overlap') !=0:
#     procinfo['nstacks'] = math.floor(procinfo.get('nstacks') *  (100/(100-procinfo.get('overlap')))) - 1
print('Time series overlap: ' + str(procinfo.get('overlap'))+'%')
print('No. of stacks: '+ str(procinfo.get('nstacks')))
print('--------------------')
print('\nBand averaging over target frequencies:')
#
#
#==================== Start band averaging ====================
#
#
# Band average value after calibration and averaging using parzen window
ftlist,bandavg = mtproc.bandavg(ts,procinfo)
#
#==================== Band averaging finished ====================

#====== Tipper calculation =========
[TxAll, TyAll] = tipper.tippall(bandavg)
timer_end = time.time()
print('\nElapsed time: ' + str(timer_end - timer_start)+'s')
del timer_start, timer_end
print('Finished.')
mahaWtTx, Tx_mcd_mean, Txmahal_robust = tipper.mcd(TxAll)
mahaWtTy, Ty_mcd_mean, Tymahal_robust = tipper.mcd(TyAll)
#====== Tipper calculated =========
# Coherency Threshold
cohMatrixEx = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
cohMatrixEy = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
mpdmat = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
AllcohEx = coherency.cohEx(bandavg)
AllcohEy = coherency.cohEy(bandavg)
alpha_degH,alpha_degE = mtproc.mpdvalues(bandavg)
#
#
ctflag = 0
if ctflag == 1:
    CohThre = [0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6]
    for i in range(np.shape(AllcohEx)[0]):
        for j in range(np.shape(AllcohEx)[1]):
            if AllcohEx[i,j] < CohThre[i]:
                cohMatrixEx[i,j] = 0
            else:
                cohMatrixEx[i,j] = 1
            if AllcohEy[i,j] < CohThre[i]:
                cohMatrixEy[i,j] = 0
            else:
                cohMatrixEy[i,j] = 1
#
mpdflag = 0
if mpdflag == 1:
    mpdlim = [-10,10]
    for i in range(np.shape(mpdmat)[0]):
        for j in range(np.shape(mpdmat)[1]):
            if alpha_degE[i,j] > mpdlim[0] and alpha_degE[i,j] < mpdlim[1]:
                mpdmat[i,j] = 0
            else:
                mpdmat[i,j] = 1
#            
mpdflag = 0
if mpdflag == 1:
    for i in range(np.shape(mpdmat)[0]):
        for j in range(np.shape(mpdmat)[1]):
            if j > 300 and j < 403:
                mpdmat[i,j] = 0
            else:
                mpdmat[i,j] = 1
    
    
bandavg['cohMatrixEx'] = cohMatrixEx
bandavg['cohMatrixEy'] = cohMatrixEy

bandavg['pre_sel_matEx'] = cohMatrixEx * mpdmat
bandavg['pre_sel_matEy'] = cohMatrixEy * mpdmat


####################

bandavg['mdmatrixEx'],bandavg['Zxx_mcd'],bandavg['Zxy_mcd'],bandavg['mahal_robustEx'] = mahaDist.mcd(bandavg,'Ex')
bandavg['mdmatrixEy'],bandavg['Zyx_mcd'],bandavg['Zyy_mcd'],bandavg['mahal_robustEy'] = mahaDist.mcd(bandavg,'Ey')
bandavg['coh_selectedEx'] = bandavg.get('mdmatrixEx')
bandavg['coh_selectedEy'] = bandavg.get('mdmatrixEy')
bandavg['tipp_selected'] = mahaWtTx * mahaWtTy
# bandavg['coh_selected'] = selmatEx * selmatEy
bandavg['avgt'] = np.sum((bandavg.get('coh_selectedEx'))!=0,axis=1)
del cohMatrixEx,cohMatrixEy
#
#
#=========== Tipper estimation ===================================
[Tx, Ty] = tipper.tipper(bandavg)
[TxVar, TyVar] = tipper.tipperVar(bandavg)
#=========== Tipper estimation DONE===============================
#
#==================== Robust estimation begins ===================
#
#
#
# #del ts
Zxx_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
Zxy_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
Zyx_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
Zyy_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
# Tx_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
# Ty_jackk = np.empty((np.shape(ftlist)[0],1),dtype=complex)
for stacki in range(np.shape(ftlist)[0]):
    print('\nComputing Jackknife estimate....ft = ' + str(ftlist[stacki,0]))
    bandavg_singleEx = mtproc.makeband(bandavg,stacki,'coh_selectedEx')
    bandavg_singleEy = mtproc.makeband(bandavg,stacki,'coh_selectedEy')
    # bandavgT_single = tipper.makeband(bandavg,stacki)
    Zxx_jackk[stacki,0],Zxy_jackk[stacki,0] = (mtproc.getjackknife(bandavg_singleEx,'Ex'))
    Zyx_jackk[stacki,0],Zyy_jackk[stacki,0] = (mtproc.getjackknife(bandavg_singleEy,'Ey'))
    # Tx_jackk[stacki,0],Ty_jackk[stacki,0] = (tipper.getjackknife(bandavgT_single))
print('Finished.')
Z_jackk = {'Zxx': Zxy_jackk}
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
    bandavg_singleEx = mtproc.makeband(bandavg,stacki,'coh_selectedEx')
    bandavg_singleEy = mtproc.makeband(bandavg,stacki,'coh_selectedEy')
    Zxx_huber[stacki,0],Zxy_huber[stacki,0],bandavgEx_huber = mtproc.huberEx(bandavg_singleEx,Z_jackk,stacki)
    Zyy_huber[stacki,0],Zyx_huber[stacki,0],bandavgEy_huber = mtproc.huberEy(bandavg_singleEy,Z_jackk,stacki)
print('Finished.')
Z_huber['Zxx'] = Zxx_huber
Z_huber['Zxy'] = Zxy_huber
Z_huber['Zyy'] = Zyy_huber
Z_huber['Zyx'] = Zyx_huber
#Z_tukey = {}
#print('\nComputing Tukey estimate....')
# Z_tukey['Zxx'],Z_tukey['Zxy'],tukey_matrixEx = mtproc.tukeyEx(bandavgEx_huber)
# Z_tukey['Zyy'],Z_tukey['Zyx'],tukey_matrixEy = mtproc.tukeyEy(bandavgEy_huber)
#print('Finished.')
Zvar = {}
Zvar['xx'],Zvar['xy'],cohEx = var.ZExvar(Z_huber,bandavg)
Zvar['yx'],Zvar['yy'],cohEy = var.ZEyvar(Z_huber,bandavg)
#
#
# ----------------------------------------------------------------------
# Apparant resistivities and phase
rho_xy = (0.2/ftlist) * ((abs(Z_huber.get('Zxy'))) ** 2)
rho_yx = (0.2/ftlist) * ((abs(Z_huber.get('Zyx'))) ** 2)
phase_xy = np.degrees(np.arctan(Z_huber.get('Zxy').imag/Z_huber.get('Zxy').real))
phase_yx = np.degrees(np.arctan(Z_huber.get('Zyx').imag/Z_huber.get('Zyx').real))
#
#
# Errors for app. resistivity and phase
err_rxy = (0.4/ftlist) * ((abs(Z_huber.get('Zxy')))) * Zvar.get('xy')
err_ryx = (0.4/ftlist) * ((abs(Z_huber.get('Zyx')))) * Zvar.get('yx')
err_pxy = np.degrees((Zvar.get('xy')) / abs(Z_huber.get('Zxy')))
err_pyx = np.degrees((Zvar.get('yx')) / abs(Z_huber.get('Zyx')))
err_rxy = err_rxy.reshape((-1,))
err_ryx = err_ryx.reshape((-1,))
err_pxy = err_pxy.reshape((-1,))
err_pyx = err_pyx.reshape((-1,))
#
# Plotting section
#
plt.figure(num=1)
plt.subplot(211)
plt.scatter(ftlist,rho_xy,c='r',s=10)
plt.scatter(ftlist,rho_yx,c='b',s=10)
plt.errorbar(ftlist,rho_xy,err_rxy,ecolor='r',fmt="none")
plt.errorbar(ftlist,rho_yx,err_ryx,ecolor='b',fmt="none")
plt.xscale('log')
plt.yscale('log')
plt.xlim((10000, 0.001)) 
plt.ylim(0.1, 100000)
plt.yticks([0.1, 1, 10, 100, 1000, 10000])
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
plt.subplot(212)
plt.scatter(ftlist,phase_xy,c='r',s=10)
plt.scatter(ftlist,phase_yx,c='b',s=10)
plt.errorbar(ftlist,phase_xy,err_pxy,ecolor='r',fmt="none")
plt.errorbar(ftlist,phase_yx,err_pyx,ecolor='b',fmt="none")
plt.xscale('log')
plt.xlim((10000, 0.001))
plt.ylim((0, 90))
plt.yticks([0,15,30,45,60,75,90])
plt.grid(which='both',linestyle='-.', linewidth=0.4)
# plt.savefig('C:/Users/ajithabh/Desktop/MTProc/'+ procinfo.get('meas')+'.png', format='png', dpi=600)
plt.figure(num=2)
plt.scatter(ftlist,cohEx,c='r')
plt.scatter(ftlist,cohEy,c='b')
plt.xscale('log')
plt.xlim((10000, 0.001)) 
plt.ylim(0, 1)
plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
# plt.savefig('C:/Users/ajithabh/Desktop/MTProc/'+ 'coh' + procinfo.get('meas')+'.png', format='png', dpi=600)
# Plot tipper =====
TxA = np.sqrt((np.real(Tx) ** 2) + (np.imag(Tx) ** 2))
TyA = np.sqrt((np.real(Ty) ** 2) + (np.imag(Ty) ** 2))
TxP = np.degrees(np.arctan2(Tx.imag,Tx.real))
TyP = np.degrees(np.arctan2(Ty.imag,Ty.real))
plt.figure(num=3)
plt.subplot(211)
plt.scatter(ftlist,TxA,c='r')
plt.scatter(ftlist,TyA,c='b')
plt.ylim(0, 1)
plt.xscale('log')
plt.xlim((10000, 0.001))
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
plt.subplot(212)
plt.scatter(ftlist,TxP,c='r')
plt.scatter(ftlist,TyP,c='b')
plt.xscale('log')
plt.xlim((10000, 0.001))
plt.grid(which='both',linestyle='-.', linewidth=0.4)
#----------
plt.figure(num=4)
plt.subplot(211)
plt.scatter(ftlist,Tx.real,c='r')
plt.scatter(ftlist,Ty.real,c='b')
#plt.ylim(0, 1)
plt.xscale('log')
plt.xlim((10000, 0.001))
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
plt.subplot(212)
plt.scatter(ftlist,Tx.imag,c='r')
plt.scatter(ftlist,Ty.imag,c='b')
plt.xscale('log')
plt.xlim((10000, 0.001))
plt.grid(which='both',linestyle='-.', linewidth=0.4)