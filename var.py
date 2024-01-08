# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:37:20 2020

@author: ajithabh
Dr. Manoj Nair helped with equations of variance and fisher distribution.

This module is used to compute variances.
The function ZExvar will compute variance for Ex component
and
The function ZEyvar will compute variance for Ey component

Both the functions will together return variances for
components Zxx,Zxy,Zyx,Zyy

The 'fisher' function is used to calculate Fisher distribution.
"""

import numpy as np
def ZExvar(Z_huber,bandavg):
    """

    Parameters
    ----------
    Z_huber : It is a Python dictionary containing huber estimates.
    bandavg : It is a Python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    ZxxVar : It is an array of float containing variance values.
    ZxyVar : It is an array of float containing variance values.
    cohEx : It is an array of float containing coherency value.

    """
    #Ex
    cohMatrixEx = bandavg.get('selectedEx')
    ExExc = bandavg.get('ExExc') * cohMatrixEx
    ExHxc = bandavg.get('ExHxc') * cohMatrixEx
    ExHyc = bandavg.get('ExHyc') * cohMatrixEx
    HxHxc = bandavg.get('HxHxc') * cohMatrixEx
    HyHxc = bandavg.get('HyHxc') * cohMatrixEx
    HxHyc = bandavg.get('HxHyc') * cohMatrixEx
    HyHyc = bandavg.get('HyHyc') * cohMatrixEx
    #
    ExExc = (np.sum(ExExc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    ExHxc = (np.sum(ExHxc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    ExHyc = (np.sum(ExHyc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    HxHxc = (np.sum(HxHxc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    HyHxc = (np.sum(HyHxc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    HxHyc = (np.sum(HxHyc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    HyHyc = (np.sum(HyHyc,axis=1)/np.sum(cohMatrixEx,axis=1)).reshape(-1,1)
    Zxx = Z_huber.get('Zxx')
    Zxy = Z_huber.get('Zxy')
    ZpZ = Zxx * np.conj(ExHxc) + Zxy * np.conj(ExHyc)
    ZpX = Zxx * HxHxc + Zxy * HyHxc
    ZpY = Zxx * HxHyc + Zxy * HyHyc
    ZpZp = Zxx * np.conj(ZpX) + Zxy * np.conj(ZpY)
    ZpZp = ZpZp * ExExc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(ZpZp[i,0])>0:
            Ccoh[i,0] = ZpZ[i,0]/np.sqrt(ZpZp[i,0])
        else:
            Ccoh[i,0] = 1+1j
    cohEx = abs(Ccoh)
    for i in range(cohEx.shape[0]):
        if cohEx[i,0] > 1.0:
            cohEx[i,0] = 1/cohEx[i,0]
    cohEx = cohEx ** 2
    fis = fisher(bandavg.get('dof'))
    d = (4/bandavg.get('dof')) * fis * (1.0 - cohEx) * ExExc
    r2xy = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(HxHxc[i,0]) * abs(HyHyc[i,0]) > 0:
            r2xy[i,0] = (HxHyc[i,0]*HyHxc[i,0]) / (HxHxc[i,0]*HyHyc[i,0])
        else:
            r2xy[i,0] = 100
    for i in range(ZpZp.shape[0]):
        if r2xy[i,0] > 0.999:
            r2xy[i,0] = 0.999
    da = np.empty((ZpZp.shape),dtype=complex)
    db = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(HxHxc[i,0]*(1-r2xy[i,0])) > 0:
            da[i,0] = d[i,0]/(HxHxc[i,0]*(1-r2xy[i,0]))
        else:
            da[i,0] = 100
    for i in range(ZpZp.shape[0]):
        if abs(HyHyc[i,0]*(1-r2xy[i,0])) > 0:
            db[i,0] = d[i,0]/(HyHyc[i,0]*(1-r2xy[i,0]))
        else:
            db[i,0] = 100
    ZxxVar = np.sqrt(np.real(da))
    ZxyVar = np.sqrt(np.real(db))
    ZxxVar = ZxxVar/2
    ZxyVar = ZxyVar/2
    cohEx = np.sqrt(cohEx)
    return ZxxVar,ZxyVar,cohEx


def ZEyvar(Z_huber,bandavg):
    """

    Parameters
    ----------
    Z_huber : It is a Python dictionary containing huber estimates.
    bandavg : It is a Python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    ZyxVar : It is an array of float containing variance values.
    ZyyVar : It is an array of float containing variance values.
    cohEy : It is an array of float containing coherency value.

    """
    cohMatrixEy = bandavg.get('selectedEy')
    EyEyc = bandavg.get('EyEyc') * cohMatrixEy
    EyHxc = bandavg.get('EyHxc') * cohMatrixEy
    EyHyc = bandavg.get('EyHyc') * cohMatrixEy
    HxHxc = bandavg.get('HxHxc') * cohMatrixEy
    HyHxc = bandavg.get('HyHxc') * cohMatrixEy
    HxHyc = bandavg.get('HxHyc') * cohMatrixEy
    HyHyc = bandavg.get('HyHyc') * cohMatrixEy
    #
    EyEyc = (np.sum(EyEyc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    EyHxc = (np.sum(EyHxc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    EyHyc = (np.sum(EyHyc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    HxHxc = (np.sum(HxHxc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    HyHxc = (np.sum(HyHxc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    HxHyc = (np.sum(HxHyc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    HyHyc = (np.sum(HyHyc,axis=1)/np.sum(cohMatrixEy,axis=1)).reshape(-1,1)
    Zyy = Z_huber.get('Zyy')
    Zyx = Z_huber.get('Zyx')
    ZpZ = Zyx * np.conj(EyHxc) + Zyy * np.conj(EyHyc)
    ZpX = Zyx * HxHxc + Zyy * HyHxc
    ZpY = Zyx * HxHyc + Zyy * HyHyc
    ZpZp = Zyx * np.conj(ZpX) + Zyy * np.conj(ZpY)
    ZpZp = ZpZp * EyEyc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(ZpZp[i,0])>0:
            Ccoh[i,0] = ZpZ[i,0]/np.sqrt(ZpZp[i,0])
        else:
            Ccoh[i,0] = 1+1j
    cohEy = abs(Ccoh)
    for i in range(cohEy.shape[0]):
        if cohEy[i,0] > 1.0:
            cohEy[i,0] = 1/cohEy[i,0]
    cohEy = cohEy ** 2
    fis = fisher(bandavg.get('dof'))
    d = (4/bandavg.get('dof')) * fis * (1.0 - cohEy) * EyEyc
    r2xy = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(HxHxc[i,0]) * abs(HyHyc[i,0]) > 0:
            r2xy[i,0] = (HxHyc[i,0]*HyHxc[i,0]) / (HxHxc[i,0]*HyHyc[i,0])
        else:
            r2xy[i,0] = 100
    for i in range(ZpZp.shape[0]):
        if r2xy[i,0] > 0.999:
            r2xy[i,0] = 0.999
    da = np.empty((ZpZp.shape),dtype=complex)
    db = np.empty((ZpZp.shape),dtype=complex)
    for i in range(ZpZp.shape[0]):
        if abs(HxHxc[i,0]*(1-r2xy[i,0])) > 0:
            da[i,0] = d[i,0]/(HxHxc[i,0]*(1-r2xy[i,0]))
        else:
            da[i,0] = 100
    for i in range(ZpZp.shape[0]):
        if abs(HyHyc[i,0]*(1-r2xy[i,0])) > 0:
            db[i,0] = d[i,0]/(HyHyc[i,0]*(1-r2xy[i,0]))
        else:
            db[i,0] = 100
    ZyxVar = np.sqrt(np.real(da))
    ZyyVar = np.sqrt(np.real(db))
    ZyxVar = ZyxVar/2
    ZyyVar = ZyyVar/2
    cohEy = np.sqrt(cohEy)
    return ZyxVar,ZyyVar,cohEy

def fisher(nue):
    """

    Parameters
    ----------
    nue : It is an array of float containing degree of freedom values.

    Returns
    -------
    w : It is an array of float containing fisher distribution values.

    """
    w = np.empty(nue.shape,dtype=float)
    for i in range(nue.shape[0]):
        if nue[i] < 6.0:
            w[i] = 5 ** (nue[i]+1.0)
        else:
            ya = 1.6446 # 5% error
            h = 2.0/(1.0/(nue[i]-5.0)+1.0/3.0)
            al = (ya * ya) / 6.0 - 0.5
            w_t = (ya * np.sqrt(h+al)/h - ((1.0/3.0) - 1.0/(nue[i]-5.0)) * 
                  (al+(5.0/6.0) - (2.0/3.0)/h))
            w[i] = np.exp(2*w_t)
    return w
