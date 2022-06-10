# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:52:42 2020

@author: AJITHABH
"""

import matplotlib
#matplotlib.use('TkAgg')
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

Z_all = bandavg.get('Zxy_single')
#coh_selected_all = bandavg.get('coh_selectedEx')
#coh_selected_all = bandavg.get('cohMatrixEx')
coh_selected_all = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
# Z_huber = Z_huber.get('Zyx')
for fnum in range(np.size(ftlist)):
    Z = Z_all[fnum,:]
    coh_selected = coh_selected_all[fnum,:].reshape(-1,1)
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
    # plt.scatter(Z_huber[fnum].real,Z_huber[fnum].imag,c='black',marker=(5, 1))
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+
              ' ('+str(procinfo.get('fs'))+' Hz) f='+ str(round(ftlist[fnum][0],2)) +' Hz')
    # plt.savefig('C:/Users/Ajithabh/Desktop/myImagePDF.eps', format='eps', dpi=1200)