# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:38:27 2020

@author: ajithabh

This script can be used to plot polarization directions
for all the time windows for all target frequencies.

X axis: Time window/Event
Y axis (Subplot 1): Magnetic polarization direction (alphaH)
Y axis (Subplot 2): Electric polarization direction (alphaE)

Number of figures generated will be equal to the number
of target frequencies. First figure will represent values
of first target frequency and so on.

See variable 'ftlist' in variable explorer to see all target
frequencies.
"""

import matplotlib
from matplotlib import pyplot as plt
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.0, 0.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
        'green':((0.0, 0.0, 0.0),
                 (0.1, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

#Magnetic field
afont = {'fontname':'Arial'}

HxHxc = bandavg.get('HxHxc')
HyHyc = bandavg.get('HyHyc')
HxHyc = bandavg.get('HxHyc')
nstacks = np.shape(HxHyc)[1]
alphaH = np.arctan(2*np.real(HxHyc)/(HxHxc-HyHyc))
alpha_degH = np.degrees(np.real(alphaH))
# MPD calculated



#Electric field

ExExc = bandavg.get('ExExc')
EyEyc = bandavg.get('EyEyc')
ExEyc = bandavg.get('ExEyc')
alphaE = np.arctan(2*np.real(ExEyc)/(ExExc-EyEyc))
alpha_degE = np.degrees(np.real(alphaE))
# MPD calculated

b = np.arange(0,nstacks,1)
yticks = np.arange(-90,100,20)

for fnum in range(np.size(ftlist)):
    directionH = alpha_degH[fnum,:]
    directionH = directionH.reshape(-1,1)
    directionE = alpha_degE[fnum,:]
    directionE = directionE.reshape(-1,1)
    #Event counter    
    b = b.reshape(-1,1)
    plt.figure(num=fnum)
    plt.subplot(211)
    plt.scatter(b,directionH)
    plt.ylim(-90, 90)
    plt.ylabel(r'$\alpha_H$ (deg.)')
    plt.xticks(**afont,fontsize=12)
    plt.yticks(yticks,**afont,fontsize=12)
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+
              ' ('+str(procinfo.get('fs'))+' Hz) f='+ str(round(ftlist[fnum][0],2)) +' Hz')
    plt.subplot(212)
    plt.scatter(b,directionE)
    plt.ylim(-90, 90)
    plt.ylabel(r'$\alpha_E$ (deg.)')
    plt.xlabel('Time window')
    plt.xticks(**afont,fontsize=12)
    plt.yticks(yticks,**afont,fontsize=12)
