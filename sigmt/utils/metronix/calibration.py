"""
Class for the metronix coil calibration
"""
import numpy as np
import scipy


class MetronixCalibration:
    """
    Class to perform Metronix specific calibration.
    """

    def __init__(
            self,
            xfft: np.ndarray,
            fft_freqs: np.ndarray,
            sensor_type: str,
            chopper_status: str,
            calibration_data: np.ndarray,
    ):
        """
        Constructor

        :param xfft: FFT values for the time series shape: (fft_length, nof_windows)
        :type xfft: np.ndarray
        :param fft_freqs: Frequencies for FFT values, 1D array.
        :type fft_freqs: np.ndarray
        :param sensor_type: Type of sensor, Eg: MFS06e
        :type sensor_type: str
        :param chopper_status: Chopper status. Either 'chopper_on' or 'chopper_off'
        :type chopper_status: str
        :param calibration_data: Calibration data as array. shape: (length(43), 3).
        :type calibration_data: np.ndarray
        """
        self.calibrated_data = None
        self.xfft = xfft
        self.fft_freqs = fft_freqs
        self.sensor_type = sensor_type
        self.chopper_status = chopper_status
        self.calibration_data = calibration_data
        self.perform_calibration()

    def perform_calibration(self):
        """
        Streamline it based on the sensor type.
        """
        if self.sensor_type == 'MFS06' or 'MFS06e':
            self.mfs06e()
        elif self.sensor_type == 'MFS07e':
            self.mfs07e()

    def mfs06e(self):
        """
        Streamline calibration for mfs06e
        """
        magnitude = self.calibration_data[:, 0] * self.calibration_data[:, 1]
        phase = np.radians(self.calibration_data[:, 2])
        calt = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
        # If Chopper is Off
        if self.chopper_status == 'chopper_off':
            interp_func = scipy.interpolate.interp1d(
                self.calibration_data[:, 0],
                calt,
                kind='linear',
                fill_value='extrapolate'
            )
            cal_all_band = interp_func(self.fft_freqs)
            self.calibrated_data = self.xfft / cal_all_band[:, np.newaxis]
        # If Chopper is On
        if self.chopper_status == 'chopper_on':
            minfindx = np.where(self.fft_freqs < 0.1)[0]
            interp_func = scipy.interpolate.interp1d(
                self.calibration_data[:, 0],
                calt,
                kind='linear',
                fill_value='extrapolate'
            )
            cal_all_band = interp_func(
                self.fft_freqs[np.max(minfindx) + 1:np.shape(self.fft_freqs)[0]])
            thmag = np.zeros(np.shape(minfindx), )
            thmag[:, ] = 0.2 * self.fft_freqs[0:np.max(minfindx) + 1]
            thph = np.arctan2(4.0, self.fft_freqs[0:np.max(minfindx) + 1])
            th_band = (thmag * np.cos(thph) + (1j * thmag * np.sin(thph))) * 1000
            cal_all_band = np.concatenate((th_band, cal_all_band))
            cal_all_band[0] = 1 + 1j  # To avoid division of zero
            self.calibrated_data = self.xfft / cal_all_band[:, np.newaxis]

    def mfs07e(self):
        """
        Streamline calibration for mfs07e
        """
        magnitude = self.calibration_data[:, 0] * self.calibration_data[:, 1]
        phase = np.radians(self.calibration_data[:, 2])
        calt = (magnitude * np.cos(phase) + (1j * magnitude * np.sin(phase))) * 1000
        # If Chopper is Off
        if self.chopper_status == 'chopper_off':
            interp_func = scipy.interpolate.interp1d(
                self.calibration_data[:, 0],
                calt,
                kind='linear',
                fill_value='extrapolate'
            )
            cal_all_band = interp_func(self.fft_freqs)
            self.calibrated_data = self.xfft / cal_all_band[:, np.newaxis]
        # If Chopper is On
        if self.chopper_status == 'chopper_on':
            minfindx = np.where(self.fft_freqs < 0.4)[0]
            interp_func = scipy.interpolate.interp1d(
                self.calibration_data[:, 0],
                calt,
                kind='linear',
                fill_value='extrapolate'
            )
            cal_all_band = interp_func(
                self.fft_freqs[np.max(minfindx) + 1:np.shape(self.fft_freqs)[0]])
            thmag = np.zeros(np.shape(minfindx), )
            thmag[:, ] = 0.2 * self.fft_freqs[0:np.max(minfindx) + 1]
            thph = np.arctan2(32.0, self.fft_freqs[0:np.max(minfindx) + 1])
            th_band = (thmag * np.cos(thph) + (1j * thmag * np.sin(thph))) * 1000
            cal_all_band = np.concatenate((th_band, cal_all_band))
            cal_all_band[0] = 1 + 1j  # To avoid division of zero
            self.calibrated_data = self.xfft / cal_all_band[:, np.newaxis]
