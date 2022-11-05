# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:41:03 2022
@author: AJITHABH K. S.
Last modified: 21-07-2022

This script can be used to plot apparent resistivity and phase value
from the text file created after processing of the data.
"""
#
import pandas as pd
import numpy as np
# from matplotlib import rc
from matplotlib import pyplot as plt
import matplotlib as mpl
# mpl.rcParams['pdf.fonttype'] = 3
path = 'D:/NGRI/PROCESSING FOLDER/AMT/TK4PON/bands/all.txt'
data = pd.read_csv(path,delimiter = "\t")
data = data.to_numpy()
ftlist = data[:,0]
#
#
Zxy = data[:,3] + 1j * data[:,4]
Zyx = data[:,5] + 1j * data[:,6]
zxyvar = data[:,10]
zyxvar = data[:,11]
# Apparant resistivities and phase
rho_xy = (0.2/ftlist) * ((abs(Zxy) ** 2))
rho_yx = (0.2/ftlist) * ((abs(Zyx) ** 2))
phase_xy = np.degrees(np.arctan(Zxy.imag/Zxy.real))
phase_yx = np.degrees(np.arctan(Zyx.imag/Zyx.real))
#
#
# Errors for app. resistivity and phase
err_rxy = (0.4/ftlist) * ((abs(Zxy))) * zxyvar
err_ryx = (0.4/ftlist) * ((abs(Zyx))) * zyxvar
err_pxy = np.degrees(zxyvar / abs(Zxy))
err_pyx = np.degrees(zyxvar / abs(Zyx))
err_rxy = err_rxy.reshape((-1,))
err_ryx = err_ryx.reshape((-1,))
err_pxy = err_pxy.reshape((-1,))
err_pyx = err_pyx.reshape((-1,))
#
afont = {'fontname':'Arial'}
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
# plt.ylim(0.01, 10000)
plt.xticks(**afont,fontsize=12)
plt.yticks(**afont,fontsize=12)
plt.xlabel('Frequency (Hz)', **afont, fontsize=12)
plt.ylabel('App. Res (Ohm.m.)',**afont, fontsize=12)
plt.yticks([0.1, 1, 10, 100, 1000, 10000])
# plt.yticks([10, 100, 1000, 10000, 100000, 1000000])
plt.xticks([10000,1000,100,10,1,0.1,0.01,0.001])
ax = plt.gca()
ax.set_box_aspect(0.5)
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.title('TK4 POWER ON')
plt.subplot(212)
plt.scatter(ftlist,phase_xy,c='r',s=10)
plt.scatter(ftlist,phase_yx,c='b',s=10)
plt.errorbar(ftlist,phase_xy,err_pxy,ecolor='r',fmt="none")
plt.errorbar(ftlist,phase_yx,err_pyx,ecolor='b',fmt="none")
plt.xscale('log')
plt.xlim((10000, 0.001)) 
plt.ylim((0, 90))
plt.xticks(**afont,fontsize=12)
plt.yticks(**afont,fontsize=12)
plt.xticks([10000,1000,100,10,1,0.1,0.01,0.001])
plt.yticks([0,15,30,45,60,75,90])
plt.xlabel('Frequency (Hz)',**afont, fontsize=12)
plt.ylabel('Phase (Degrees)',**afont, fontsize=12)
ax = plt.gca()
ax.set_box_aspect(0.5)
plt.grid(which='both',linestyle='-.', linewidth=0.4)
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.6)
