"""
@author: AJITHABH K. S.
Last modified: 27-07-2022

This module contains functions to perform calculations for
different data selection methods.
"""
#######################################################################
"""
This funtions cohEx and cohEy compute the coherency values for all 
time windows with respect to corresponding target frequencies.

The inputs are the auto- and cross-spectra values and impedance
values for all time windows with respect to corresponding target frequencies.

The outputs are coherency values for Ex and Ey components. 

The variable AllcohEx contains coherency values for Ex components. Number of rows represent
number of target frequencies and column represent number of time windows/events.

The variable AllcohEy contains coherency values for Ey components. Number of rows represent
number of target frequencies and column represent number of time windows/events.
"""
import numpy as np
def cohEx(bandavg):
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

"""
This funtions returns magnetic polarization direction (alpha_degH) and electric
polarization direction (alpha_degE) for all time windows for all target frequencies.
Number of rows represent number of target frequencies and column represent number 
of time windows/events.
"""

def pdvalues(bandavg):
    HxHxc = bandavg.get('HxHxc')
    HyHyc = bandavg.get('HyHyc')
    HxHyc = bandavg.get('HxHyc')
    nstacks = np.shape(HxHyc)[1]
    alphaH = np.arctan(2*np.real(HxHyc)/(HxHxc-HyHyc))
    alpha_degH = np.degrees(np.real(alphaH))

    #Electric field
    
    ExExc = bandavg.get('ExExc')
    EyEyc = bandavg.get('EyEyc')
    ExEyc = bandavg.get('ExEyc')
    alphaE = np.arctan(2*np.real(ExEyc)/(ExExc-EyEyc))
    alpha_degE = np.degrees(np.real(alphaE))
    # PD calculated
    
    return alpha_degH,alpha_degE
