"""
Class for the Phoenix coil calibration

"""
from typing import Dict

import numpy as np
import scipy


class PhoenixCalibration:
    """
    Class to perform Phoenix specific calibration.
    """

    def __init__(
            self,
            xfft: np.ndarray,
            fft_freqs: np.ndarray,
            calibration_data: Dict
    ):
        """
        Constructor

        :param xfft: FFT values for the time series shape: (fft_length, nof_windows)
        :type xfft: np.ndarray
        :param fft_freqs: Frequencies for FFT values, 1D array.
        :type fft_freqs: np.ndarray
        :param calibration_data: Calibration data as dictionary
        :type calibration_data: Dict
        """
        self.calibrated_data = None
        self.xfft = xfft
        self.fft_freqs = fft_freqs
        self.calibration_data = calibration_data
        self.perform_calibration()

    def perform_calibration(self):
        """
        Apply calibration
        """
        magnitude = np.array(self.calibration_data['magnitude'])
        phase = np.radians(self.calibration_data['phs_deg'])
        calibration = (magnitude * np.cos(phase)) + (1j * magnitude * np.sin(phase))
        interp_func = scipy.interpolate.interp1d(
            np.array(self.calibration_data['freq_Hz']),
            calibration,
            kind='linear',
            fill_value='extrapolate'
        )
        cal_all_band = interp_func(self.fft_freqs)
        self.calibrated_data = self.xfft / cal_all_band[:, np.newaxis]
