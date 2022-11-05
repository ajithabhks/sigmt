"""
@author: AJITHABH K. S.
Last modified: 27-07-2022

This module contains functions to perform calculations for
different data selection methods.
"""
#######################################################################
import numpy as np
def cohEx(bandavg):
    """
    
    Parameters
    ----------
    bandavg : It is a Python dictionary containing the auto- and cross- spectra
        values and impedance values for all time windows at all target frequencies. 

    Returns
    -------
    AllcohEx : It is an array containing coherency values of Ex component for
        all time windows for all target frequencies. Number of rows represent
        number of target frequencies and number of column represent number of time 
        windows/events.

    """
    #Ex predicted
    ExExc = bandavg.get('ExExc')
    ExHxc = bandavg.get('ExHxc')
    ExHyc = bandavg.get('ExHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHxc = bandavg.get('HyHxc')
    HxHyc = bandavg.get('HxHyc')
    HyHyc = bandavg.get('HyHyc')
    Zxx = bandavg.get('Zxx_single')
    Zxy = bandavg.get('Zxy_single')
    ZpZ = Zxx * np.conj(ExHxc) + Zxy * np.conj(ExHyc)
    ZpX = Zxx * HxHxc + Zxy * HyHxc
    ZpY = Zxx * HxHyc + Zxy * HyHyc
    ZpZp = Zxx * np.conj(ZpX) + Zxy * np.conj(ZpY)
    ZpZp = ZpZp * ExExc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for j in range(ZpZp.shape[0]):
        for i in range(ZpZp.shape[1]):
            if abs(ZpZp[j,i])>0:
                Ccoh[j,i] = ZpZ[j,i]/np.sqrt(ZpZp[j,i])
            else:
                Ccoh[j,i] = 1+1j
        cohEx = abs(Ccoh)
        for i in range(cohEx.shape[1]):
            if cohEx[j,i] > 1.0:
                cohEx[j,i] = 1/cohEx[j,i]
    AllcohEx = cohEx
    return AllcohEx


def cohEy(bandavg):
    """
    
    Parameters
    ----------
    bandavg : It is a Python dictionary containing the auto- and cross- spectra
        values and impedance values for all time windows at all target frequencies. 

    Returns
    -------
    AllcohEx : It is an array containing coherency values of Ey component for
        all time windows for all target frequencies. Number of rows represent
        number of target frequencies and number of column represent number of time 
        windows/events.

    """
    #Ey predicted
    EyEyc = bandavg.get('EyEyc')
    EyHxc = bandavg.get('EyHxc')
    EyHyc = bandavg.get('EyHyc')
    HxHxc = bandavg.get('HxHxc')
    HyHxc = bandavg.get('HyHxc')
    HxHyc = bandavg.get('HxHyc')
    HyHyc = bandavg.get('HyHyc')
    Zyy = bandavg.get('Zyy_single')
    Zyx = bandavg.get('Zyx_single')
    ZpZ = Zyx * np.conj(EyHxc) + Zyy * np.conj(EyHyc)
    ZpX = Zyx * HxHxc + Zyy * HyHxc
    ZpY = Zyx * HxHyc + Zyy * HyHyc
    ZpZp = Zyx * np.conj(ZpX) + Zyy * np.conj(ZpY)
    ZpZp = ZpZp * EyEyc
    Ccoh = np.empty((ZpZp.shape),dtype=complex)
    for j in range(ZpZp.shape[0]):
        for i in range(ZpZp.shape[1]):
            if abs(ZpZp[j,i])>0:
                Ccoh[j,i] = ZpZ[j,i]/np.sqrt(ZpZp[j,i])
            else:
                Ccoh[j,i] = 1+1j
        cohEy = abs(Ccoh)
        #cohEy = 1-cohEy
        for i in range(cohEy.shape[1]):
            if cohEy[j,i] > 1.0:
                cohEy[j,i] = 1/cohEy[j,i]
    AllcohEy = cohEy
    return AllcohEy


def pdvalues(bandavg):
    """
    
    Parameters
    ----------
    bandavg : It is a Python dictionary containing the auto- and cross- spectra
        values and impedance values for all time windows at all target frequencies.

    Returns
    -------
    alpha_degH : It is an array containing magnetic polarization directions for 
        all time windows at all target frequencies. Number of rows represent number 
        of target frequencies and number of column represent number of time windows/events.
    alpha_degE : It is an array containing electric polarization direction for 
        all time windows at all target frequencies. Number of rows represent number 
        of target frequencies and number of column represent number of time windows/events.

    """
    HxHxc = bandavg.get('HxHxc')
    HyHyc = bandavg.get('HyHyc')
    HxHyc = bandavg.get('HxHyc')
    nstacks = np.shape(HxHyc)[1]
    alphaH = np.arctan(2*np.real(HxHyc)/(HxHxc-HyHyc))
    alpha_degH = np.degrees(np.real(alphaH))

    ExExc = bandavg.get('ExExc')
    EyEyc = bandavg.get('EyEyc')
    ExEyc = bandavg.get('ExEyc')
    alphaE = np.arctan(2*np.real(ExEyc)/(ExExc-EyEyc))
    alpha_degE = np.degrees(np.real(alphaE))
    return alpha_degH,alpha_degE
