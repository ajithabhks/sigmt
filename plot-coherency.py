# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:52:42 2020

@author: ajithabh

This script can be used to plot coherency values
for all the time windows for all target frequencies.

This script enables ploting of the coherency values
for Ex and Ey components in an argand diagram. Any impedance
value can be selected for plotting. The color of the data points
indicate the coherency value.

To plot Ex component, make 
cc = AllcohEx[fnum,:]

Similarly for Ey, make
cc = AllcohEy[fnum,:]

Change 'Zxy_single' to 'Zyx_single' in following line
to choose argand diagram of Zyx component.
Z_all = bandavg.get('Zxy_single')

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
"""
Try 'Zxy_single', 'Zyx_single', 'Zxx_single', 'Zyy_single'
in the below field as required
"""
Z_all = bandavg.get('Zxy_single')
coh_selected_all = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
for fnum in range(np.size(ftlist)):
    Z = Z_all[fnum,:]
    coh_selected = coh_selected_all[fnum,:].reshape(-1,1)
    """
    Try: cc = AllcohEx[fnum,:] to plot Ex component
    Try: cc = AllcohEy[fnum,:] to plot Ey component
    """
    cc = AllcohEx[fnum,:]
    cc = cc.reshape(-1,1)
    Z = Z.reshape(-1,1)
    ind_coh = np.where(coh_selected==0)[0].reshape(-1,1)
    cc = np.delete(cc,ind_coh).reshape(-1,1)
    Z = np.delete(Z,ind_coh).reshape(-1,1)
    Z_real = np.real(Z)
    Z_imag = np.imag(Z)
    plt.figure(num=fnum)
    sc = plt.scatter(Z_real,Z_imag,c=cc,cmap=my_cmap)
    plt.colorbar(sc)
    plt.clim(0,1) 
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz) f='+ str(round(ftlist[fnum][0],2)) +' Hz')
