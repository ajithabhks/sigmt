# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 12:41:28 2021
@author: AJITHABH K. S.
Last modified: 13-07-2022

This module contains functions to calculate Mahalanobis distance (MD) for every
time windows for all target frequencies. It is also used to select time windows/events within
a MD threshold.

Input files are impedance values and a matrix contains pre-selection information.
The discarded time windows will have value '0' and selected windows will have value
'1' in the pre-selection matrix.

MD will be calculated for every selected time windows/events.
An MD value of 10 will be assigned for the discarded time windows.

The output of the module consists of a matrix contains selected events details (mahaWt).
The value of '1' is given for the events falls within the MD threshold and '0' for others.
The MD thresholds can be set in the 'config.py' file.

Other outputs are:-
Z1_mcd_mean: Robust cluster center
Z2_mcd_mean: Robust cluster center
mahal_robust: MD values for all events

"""

def mcd(bandavg,mode,config):
    """
    
    Parameters
    ----------
    bandavg : It is a python dictionary containing the auto- and cross- spectra
        values, impedance values, arrays containing pre-selection information 
        (pre_sel_matEx and pre_sel_matEy) for all time windows at all target 
        frequencies. The discarded time windows will have value '0' and selected 
        windows will have value '1' in the pre-selection arrays.
    mode : It is a string which is either 'Ex' or 'Ey'.
    config : It is a python dictionary containing processing parameters such as
        FFT length, Parzen window radius, Mahalanobis distance threshold values.

    Returns
    -------
    mahaWt : It is an array containing selected events details (mahaWt). The value 
        of '1' is given for the events falls within the MD threshold and '0' for others.
    Z1_mcd_mean : It is an array containing robust cluster centers.
    Z2_mcd_mean : It is an array containing robust cluster centers.
    mahal_robust : It is an array containing mahalanobis distances for all time 
    windows at all target frequencies.

    """
    import numpy as np
    if (mode == 'Ex'):
        Z1 = bandavg.get('Zxx_single') * bandavg.get('pre_sel_matEx')
        Z2 = bandavg.get('Zxy_single') * bandavg.get('pre_sel_matEx')
        # coh_select = bandavg.get('coh_selectedEx')
    elif (mode == 'Ey'):
        Z1 = bandavg.get('Zyx_single') * bandavg.get('pre_sel_matEy')
        Z2 = bandavg.get('Zyy_single') * bandavg.get('pre_sel_matEy')
    mahaWt = np.ones(Z1.shape)
    mahal_robust = np.empty((Z1.shape),dtype=float)
    Z1_mcd_mean = np.empty((Z1.shape[0],1),dtype=complex)
    Z2_mcd_mean = np.empty((Z1.shape[0],1),dtype=complex)
    Z1 = np.transpose(Z1)
    Z2 = np.transpose(Z2)
    for i in range(Z1.shape[1]):
        Z1_single = Z1[:,i].reshape(-1,1)
        Z2_single = Z2[:,i].reshape(-1,1)
        [mahaWt_single,Z1_mean,Z2_mean,mahal_single] = getmahaWt(Z1_single,Z2_single,config)
        Z1_mcd_mean[i] = Z1_mean
        Z2_mcd_mean[i] = Z2_mean
        for k in range(mahaWt_single.shape[0]):
            mahaWt[i,k] = mahaWt_single[k]
            mahal_robust[i,k] = mahal_single[k]
    return mahaWt, Z1_mcd_mean, Z2_mcd_mean, mahal_robust
    #
    
def getmahaWt(Z1_single,Z2_single,config):
    """
    
    Parameters
    ----------
    Z1_single : It is an array containing impedance value (either Zxx or Zyx) for 
    all time windows at a particular target frequency.
    Z2_single : It is an array containing impedance value (either Zxy or Zyy) for 
    all time windows at a particular target frequency.
    config : It is a python dictionary containing processing parameters such as
        FFT length, Parzen window radius, Mahalanobis distance threshold values.

    Returns
    -------
    mahaWt : It is an array containing selected events details (mahaWt). The value 
        of '1' is given for the events falls within the MD threshold and '0' for others.
    Z1_mean : It is an array containing robust cluster center.
    Z2_mean : It is an array containing robust cluster center.
    mahaWt_temp : It is an array containing mahalanobis distances for all time 
    windows at a particular target frequency.

    """
    import numpy as np
    nozeromat = np.where(Z1_single != 0)[0]
    zeromat = np.where(Z1_single == 0)[0]
    mahaWt_temp = np.zeros(Z1_single.shape)
    Z1_single = Z1_single[Z1_single != 0]
    Z2_single = Z2_single[Z2_single != 0]
    Z1_singleR = Z1_single.real
    Z1_singleI = Z1_single.imag
    Z2_singleR = Z2_single.real
    Z2_singleI = Z2_single.imag
    X = np.vstack((Z1_singleR, Z1_singleI,Z2_singleR, Z2_singleI))
    X = np.transpose(X)
    from sklearn.covariance import EmpiricalCovariance, MinCovDet
    # fit a MCD robust estimator to data
    robust_cov = MinCovDet().fit(X)
    mahal_robust = np.sqrt(robust_cov.mahalanobis(X))
    for m in range(nozeromat.shape[0]):
        mahaWt_temp[nozeromat[m]] = mahal_robust[m]
    for n in range(zeromat.shape[0]):
        mahaWt_temp[zeromat[n]] = 10
    mahaWt = np.zeros(([nozeromat.shape[0] + zeromat.shape[0],1]),dtype=int)
    for k in range(nozeromat.shape[0] + zeromat.shape[0]):
        if (mahaWt_temp[k] <= config.get('MD_threshold_impedance')): # MD threshold
            mahaWt[k] = 1
    Z1_mean = np.complex(robust_cov.location_[0],robust_cov.location_[1])
    Z2_mean = np.complex(robust_cov.location_[2],robust_cov.location_[3])
    return mahaWt, Z1_mean, Z2_mean, mahaWt_temp
