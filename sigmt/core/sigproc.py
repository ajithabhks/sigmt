"""
Module for signal processing operations
"""

import numpy as np
from scipy import signal


def notchfilsos(ts: np.ndarray, fs: float, notchfreq: float) -> np.ndarray:
    """
    SOS filter for notch filtering. Given frequency and all harmonics.

    :param ts: Time series data
    :type ts: np.ndarray
    :param fs: Sampling frequency
    :type fs: float
    :param notchfreq: Frequency to be filtered
    :type notchfreq: float

    :return: Filtered time series data
    :rtype: np.ndarray
    """
    min_fs = (notchfreq + 5) * 2
    if fs > min_fs:
        har = int((fs / 2.5) / notchfreq) + 1
        if har > int(14000 / notchfreq):
            har = int(14000 / notchfreq)
        print('No. of harmonics: ' + str(har))
        sos = signal.butter(4, [notchfreq - 5, notchfreq + 5], btype='bandstop', fs=fs, output='sos')
        for n in range(2, har):
            f0 = n * notchfreq
            sos_new = signal.butter(4, [f0 - 5, f0 + 5], btype='bandstop', fs=fs, output='sos')
            sos = np.concatenate((sos, sos_new), axis=0)
        ts = signal.sosfilt(sos, ts, axis=0)
    else:
        print("Skipping notch filter as it cannot be performed for this sampling frequency")

    return ts


def do_fft(ts: np.ndarray, fs: float, fft_length: int) -> tuple:
    """
    Function to perform FFT.

    :param ts: Time series data
    :type ts: np.ndarray
    :param fs: Sampling frequency
    :type fs: float
    :param fft_length: FFT length
    :type fft_length: int

    :return: Tuple of f[np.ndarray: shape n,1], xfft[np.ndarray: shape n,n]
    :rtype: tuple
    """
    w = np.hanning(fft_length).reshape(-1, 1)
    fft_value = np.fft.fft(ts * w, fft_length, axis=0)
    xfft = np.asarray(fft_value[0:int(fft_length / 2), :])
    fline = np.asarray([np.linspace(0, int(fft_length / 2), num=int(fft_length / 2), dtype=int)])
    f = np.asarray(fline * fs / fft_length).T
    f = f[:, 0]

    return f, xfft
