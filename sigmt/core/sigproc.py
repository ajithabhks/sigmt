"""
Module for signal processing operations
"""

from typing import Optional

import numpy as np
from scipy import signal


def notch_filter_sos(time_series: np.ndarray,
                     sampling_frequency: float,
                     notch_frequency: float,
                     harmonics: Optional[int] = None) -> np.ndarray:
    """
    SOS filter for notch filtering. Filters out the requested notch frequency and harmonics.

    :param time_series: Time series data
    :type time_series: np.ndarray
    :param sampling_frequency: Sampling frequency
    :type sampling_frequency: float
    :param notch_frequency: Frequency to be filtered
    :type notch_frequency: float
    :param harmonics: Number of harmonics to be filtered out.
    :type harmonics: Optional[int]

    :return: Filtered time series data
    :rtype: np.ndarray

    """
    min_fs = (notch_frequency + 5) * 2
    if sampling_frequency > min_fs:
        if harmonics is None:
            harmonics = int((sampling_frequency / 2.5) / notch_frequency) + 1
        if harmonics > int(14000 / notch_frequency):
            harmonics = int(14000 / notch_frequency)
        print('No. of harmonics: ' + str(harmonics))
        sos = signal.butter(4, [notch_frequency - 5, notch_frequency + 5], btype='bandstop',
                            fs=sampling_frequency,
                            output='sos')
        for n in range(2, harmonics + 1):
            f0 = n * notch_frequency
            high = f0 + 5
            nyq = f0/2

            # stop once beyond nyq
            if high >= nyq:
                print(f'Could not apply notch filter at frequency {f0} Hz.\n')
                break

            sos_new = signal.butter(4, [f0 - 5, f0 + 5], btype='bandstop', fs=sampling_frequency,
                                    output='sos')
            sos = np.concatenate((sos, sos_new), axis=0)
        time_series = signal.sosfilt(sos, time_series, axis=0)
    else:
        print("Skipping notch filter as it cannot be performed for this sampling frequency")

    return time_series


def do_fft(
        ts: np.ndarray,
        fs: float,
        fft_length: int
) -> tuple:
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
