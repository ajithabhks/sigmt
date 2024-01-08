# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:38:58 2020

@author: AJITHABH K. S.
Last modified: 27-07-2022

This module deals with the calculation for tipper estimation.
"""

import numpy as np
def tippall(bandavg):
    """
    This function computes tipper values for all time windows for all target 
    frequencies.

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    Tx : It is an array of complex which contains TipperX values for all time windows
        for all target frequencies.
    Ty : It is an array of complex which contains TipperY values for all time windows
        for all target frequencies.

    """
    HzHxc = bandavg.get('HzHxc')
    HzHyc = bandavg.get('HzHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHyc = bandavg.get('HyHyc')
    HxHyc = bandavg.get('HxHyc')
    HyHxc = bandavg.get('HyHxc')
    det = (HxHxc * HyHyc) - (HxHyc * HyHxc)
    Tx = ((HzHxc * HyHyc) - (HzHyc * HyHxc))/det
    Ty = ((HzHyc * HxHxc) - (HzHxc * HxHyc))/det
    return Tx,Ty


def tipper(bandavg):
    """
    This function estimates the tipper values using the selected time windows.

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    Tx : It is an array of complex which contains TipperX values for all target frequencies.
    Ty : It is an array of complex which contains TipperY values for all target frequencies.

    """
    selMatrix = bandavg.get('tipp_selected')
    HzHxc = bandavg.get('HzHxc') * selMatrix
    HzHyc = bandavg.get('HzHyc') * selMatrix
    HxHxc = bandavg.get('HxHxc') * selMatrix
    HyHyc = bandavg.get('HyHyc') * selMatrix
    HxHyc = bandavg.get('HxHyc') * selMatrix
    HyHxc = bandavg.get('HyHxc') * selMatrix
    HzHxc = (np.sum(HzHxc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    HzHyc = (np.sum(HzHyc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    HxHxc = (np.sum(HxHxc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    HyHyc = (np.sum(HyHyc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    HxHyc = (np.sum(HxHyc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    HyHxc = (np.sum(HyHxc,axis=1)/np.sum(selMatrix,axis=1)).reshape(-1,1)
    det = (HxHxc * HyHyc) - (HxHyc * HyHxc)
    Tx = 1 * ((HzHxc * HyHyc) - (HzHyc * HyHxc))/det
    Ty = 1 * ((HzHyc * HxHxc) - (HzHxc * HxHyc))/det
    return Tx,Ty

def tipperVar(bandavg):
    """
    This function computes the tipper variances.

    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.

    Returns
    -------
    TxVar : It is an array of complex which contains TipperX variance values for 
        all target frequencies.
    TyVar : It is an array of complex which contains TipperY variance values for 
        all target frequencies.

    """
    coh_selected_all = bandavg.get('tipp_selected')
    HzHxc = bandavg.get('HzHxc')
    HzHyc = bandavg.get('HzHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHyc = bandavg.get('HyHyc')
    HxHyc = bandavg.get('HxHyc')
    HyHxc = bandavg.get('HyHxc')
    det = (HxHxc * HyHyc) - (HxHyc * HyHxc)
    Tx = ((HzHxc * HyHyc) - (HzHyc * HyHxc))/det
    Ty = ((HzHyc * HxHxc) - (HzHxc * HxHyc))/det
    TxVar = np.empty((Tx.shape[0],1),dtype=float)
    TyVar = np.empty((Tx.shape[0],1),dtype=float)
    for fnum in range(np.shape(Tx)[0]):
        T1 = Tx[fnum,:]
        T2 = Ty[fnum,:]
        coh_selected = coh_selected_all[fnum,:].reshape(-1,1)
        ind_coh = np.where(coh_selected==0)[0].reshape(-1,1)
        T1 = T1.reshape(-1,1)
        T1 = np.delete(T1,ind_coh).reshape(-1,1)
        T1_std = np.std(T1)
        T2 = T2.reshape(-1,1)
        T2 = np.delete(T2,ind_coh).reshape(-1,1)
        T2_std = np.std(T2)
        TxVar[fnum] = T1_std ** 2
        TyVar[fnum] = T2_std ** 2
    return TxVar,TyVar

def mcd(T,config):
    """
    This function is used to select tipper data for final
    estimation using Mahalanobis distance.

    Parameters
    ----------
    T : It is an array of complex which contains tipper values for all time windows
        for all target frequencies.
    config : It is a python dictionary containing processing parameters such as
        FFT length, Parzen window radius, Mahalanobis distance threshold values.

    Returns
    -------
    mahaWt : It is an array of float containing zeros and ones.
    T_mcd_mean : It is an array of complex containing tipper cluster centers.

    """
    import numpy as np
    mahaWt = np.ones(T.shape)
    mahal_robust = np.empty((T.shape),dtype=float)
    T_mcd_mean = np.empty((T.shape[0],1),dtype=complex)
    T = np.transpose(T)
    for i in range(T.shape[1]):
        T_single = T[:,i].reshape(-1,1)
        [mahaWt_single,T_mean,mahal_single] = getmahaWt(T_single,config)
        T_mcd_mean[i] = T_mean
        for k in range(mahaWt_single.shape[0]):
            mahaWt[i,k] = mahaWt_single[k]
            mahal_robust[i,k] = mahal_single[k]
    return mahaWt, T_mcd_mean

def getmahaWt(T_single,config):
    """

    Parameters
    ----------
    T_single : It is an array of complex containing tipper information for a 
        target frequency.
    config : It is a python dictionary containing processing parameters such as
        FFT length, Parzen window radius, Mahalanobis distance threshold values.

    Returns
    -------
    mahaWt : It is an array of float containing zeros and ones.
    T_mean : It is an array of complex containing tipper cluster centers.
    mahal_robust : It is an array of float containing mahalanobis distances.

    """
    import numpy as np
    T_singleR = T_single.real
    T_singleI = T_single.imag
    X = np.hstack((T_singleR, T_singleI))
    from sklearn.covariance import EmpiricalCovariance, MinCovDet
    # fit a MCD robust estimator to data
    robust_cov = MinCovDet().fit(X)
    mahal_robust = np.sqrt(robust_cov.mahalanobis(X))
    mahaWt = np.zeros((mahal_robust.shape),dtype=int)
    for k in range(mahal_robust.shape[0]):
        if (mahal_robust[k] <= config.get('MD_threshold_tipper')):
            mahaWt[k] = 1.0
    #out_mahal = np.where(mahal_robust_cov > 1.5)
    T_mean = complex(robust_cov.location_[0],robust_cov.location_[1])
    return mahaWt, T_mean, mahal_robust
