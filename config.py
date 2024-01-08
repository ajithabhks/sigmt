# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 11:46:56 2022

@author: ajithabh

This is the configuration file to change FFT Length, Parzen Radius,
Mahalanobis Distance (MD) threshold for impedance and tipper calculations.
"""

def configuration():
    config = {}
    config['FFT_Length'] = 0
    # Give '0' to automatically select an FFT length according to
    # Borah et al. (2015)
    # Otherwise try 256, 512, 1024, 2048, 4096, 8192, 16384
    config['parzen_radius'] = 0
    # Give '0' to select default parzen radius
    config['MD_threshold_impedance'] = 1.5
    # Mahalanobis distance threshold value for impedance calculations
    config['MD_threshold_tipper'] = 1.0
    # Mahalanobis distance threshold value for tipper calculations
    return config