"""
Module containing statistics functions
"""

import numpy as np
import xarray as xr


def parzen(f: np.ndarray, ft: float, cr: float) -> np.ndarray:
    """
    Function to create a parzen window

    :param f: Array of frequencies from FFT
    :type f: np.ndarray
    :param ft: Target frequency
    :type ft: float
    :param cr: Parzen Window Radius
    :type cr: float
    :return: Parzen window. Length: FFT_length/2. Array shape: (n,1) - 1D
    :rtype: np.ndarray

    """
    pf = np.empty((np.shape(f)[0],))
    pf[:] = np.nan
    fr = cr * ft
    pf[0] = 0
    for i in range(1, np.shape(f)[0]):
        cond = abs(ft - f[i])
        if cond == 0:
            pf[i] = 1
        elif (cond > 0) and (cond < fr):
            u = (np.pi * abs(cond)) / fr
            pf[i] = (np.sin(u) / u) ** 4
        elif (cond > fr) or (cond == fr):
            pf[i] = 0
    pf = pf.reshape(-1, 1)
    return pf


def jackknife(z: xr.DataArray) -> complex:
    """
    Function to compute jackknife mean

    :param z: Data array of impedance values
    :type z: xr.DataArray
    :return: Jackknife mean
    :rtype: complex

    """
    n = len(z)
    jackknife_means = np.zeros(n, dtype=complex)

    for i in range(n):
        jackknife_sample = np.delete(z, i)
        jackknife_means[i] = np.mean(jackknife_sample)

    jackknife_mean = np.mean(jackknife_means)

    return jackknife_mean


def fisher_critical(dof: int) -> float:
    """
    Calculate the critical value of the F-distribution for a given
    degrees of freedom.
    Dr. Manoj C. Nair helped with this computation.

    :param dof: Degree of freedom
    :type dof: int
    :return: The critical value of F-distribution
    :rtype: float

    """

    if dof < 6.0:
        # For small degrees of freedom, use a simple exponential approximation
        f_critical = 5 ** (dof + 1.0)
    else:
        # For larger degrees of freedom, use a more accurate approximation
        z = 1.6446
        h = 2.0 / (1.0 / (dof - 5.0) + 1.0 / 3.0)
        alpha = (z * z) / 6.0 - 0.5
        w_t = (z * np.sqrt(h + alpha) / h -
               ((1.0 / 3.0) - 1.0 / (dof - 5.0)) *
               (alpha + (5.0 / 6.0) - (2.0 / 3.0) / h))
        f_critical = np.exp(2 * w_t)

    return f_critical
