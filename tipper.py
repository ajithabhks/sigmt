# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:38:58 2020

@author: AJITHABH K. S.
Last modified: 27-07-2022

This module deals with the calculation for tipper estimation.
More details of each functions are given as comments below.
"""

import numpy as np

# This function computes tipper values for all time windows
# for all target frequencies
def tippall(bandavg):
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

# This function estimates the tipper values using the selected time windows
def tipper(bandavg):
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

# This function computes the tipper variances
def tipperVar(bandavg):
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

# This was written for selective stacking for tipper estimates
# It is no longer used
def selstack(T):
    Z = T
    selmat = np.empty(Z.shape,dtype=int)
    for j in range(selmat.shape[0]):
        Z_s = Z[j,].reshape(-1,1)
        Z_s_mean = np.mean(Z_s)
        Z_s_std = np.std(Z_s)
        for i in range(selmat.shape[1]):
            if Z_s[i] >= (Z_s_mean - Z_s_std) and Z_s[i] <= (Z_s_mean + Z_s_std):
                selmat[j,i] = 1
            else:
                selmat[j,i] = 0
    return selmat

# This function is no longer used
def makeband(bandavg,i):
    HzHxc = bandavg.get('HzHxc')[i,:].reshape(1,-1)
    HzHyc = bandavg.get('HzHyc')[i,:].reshape(1,-1)
    HxHxc = bandavg.get('HxHxc')[i,:].reshape(1,-1)
    HyHyc = bandavg.get('HyHyc')[i,:].reshape(1,-1)
    HxHyc = bandavg.get('HxHyc')[i,:].reshape(1,-1)
    HyHxc = bandavg.get('HyHxc')[i,:].reshape(1,-1)
    coh_selected = bandavg.get('tipp_selected')[i,:].reshape(1,-1)
    ind_coh = np.where(coh_selected==0)[1].reshape(1,-1)
    HzHxc = np.delete(HzHxc,ind_coh).reshape(1,-1)
    HzHyc = np.delete(HzHyc,ind_coh).reshape(1,-1)
    HxHxc = np.delete(HxHxc,ind_coh).reshape(1,-1)
    HyHyc = np.delete(HyHyc,ind_coh).reshape(1,-1)
    HxHyc = np.delete(HxHyc,ind_coh).reshape(1,-1)
    HyHxc = np.delete(HyHxc,ind_coh).reshape(1,-1)
    bandavg_single = {}
    bandavg_single['HzHxc'] = HzHxc
    bandavg_single['HzHyc'] = HzHyc
    bandavg_single['HxHxc'] = HxHxc
    bandavg_single['HyHyc'] = HyHyc
    bandavg_single['HxHyc'] = HxHyc
    bandavg_single['HyHxc'] = HyHxc
    return bandavg_single

# This function is used to select tipper data for final
# estimation using Mahalanobis distance.
def mcd(T,config):
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
    #
def getmahaWt(T_single,config):
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
    T_mean = np.complex(robust_cov.location_[0],robust_cov.location_[1])
    return mahaWt, T_mean, mahal_robust
