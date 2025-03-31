"""
Class for band averaging. Here, auto and cross spectra are calculated for
each target frequencies and transfer functions are computed for all events.
"""

import gc
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

import numpy as np
import xarray as xr
from scipy import signal

import sigmt.core.sigproc as sp
import sigmt.core.statistics as stats
import sigmt.utils.metronix.calibration as calibration
import sigmt.utils.utils as utils


class BandAvg:
    """
    Class to perform band averaging.

    """

    def __init__(self,
                 time_series: dict[str, np.ndarray],
                 sampling_frequency: float,
                 fft_length: int = 1024,
                 parzen_window_radius: float = 0.25,
                 overlap: int = 50,
                 frequencies_per_decade: int = 12,
                 remote_reference: bool = False,
                 calibrate_electric: bool = False,
                 calibrate_magnetic: bool = False,
                 calibration_data_electric: Optional[dict[str, dict[str, float]]] = None,
                 calibration_data_magnetic: Optional[dict] = None,
                 apply_notch_filter: bool = False,
                 notch_frequency: Optional[float] = None,
                 process_mt: bool = True,
                 process_tipper: bool = True) -> None:
        """
        Constructor

        :param time_series: A dictionary containing time series data as 1D numpy arrays.
                    Each key represents a component of the data, and the corresponding
                    value is a numpy array of time series values. The dictionary should
                    include the following keys and associated units:
                    - 'ex': Electric field component in the x-direction
                        - Unit: mV/km (if calibrated), mV (Metronix, if not calibrated)
                    - 'ey': Electric field component in the y-direction
                        - Unit: mV/km (if calibrated), mV (Metronix, if not calibrated)
                    - 'hx': Magnetic field component in the x-direction
                        - Unit: nT (if calibrated), mV (Metronix, if not calibrated)
                    - 'hy': Magnetic field component in the y-direction
                        - Unit: nT (if calibrated), mV (Metronix, if not calibrated)
                    - 'hz': Magnetic field component in the z-direction
                        - Unit: nT (if calibrated), mV (Metronix, if not calibrated)
                    - 'rx': Remote magnetic field component in the x-direction
                        - Unit: nT (if calibrated), mV (Metronix, if not calibrated)
                    - 'ry': Remote magnetic field component in the y-direction
                        - Unit: nT (if calibrated), mV (Metronix, if not calibrated)

        :type time_series: dict[str, numpy.ndarray]
        :param sampling_frequency: Sampling frequency of the time_series.
        :type sampling_frequency: float
        :param fft_length: The length of the Fast Fourier Transform (FFT) to be applied.
                           The default value is 1024.
        :type fft_length: int
        :param parzen_window_radius: The radius of the Parzen window used for band averaging.
                                     The default value is 0.25.
        :type parzen_window_radius: float
        :param overlap: The percentage of overlap between consecutive time windows when segmenting the time series.
                        The default value is 50%.
        :type overlap: int
        :param frequencies_per_decade: The number of target frequencies to be used per decade for band averaging.
                                       The default value is 12.
        :type frequencies_per_decade: int
        :param remote_reference: A boolean flag indicating whether remote reference is used for processing.
                                 Set to True if a remote reference is used. The default is False.
        :type remote_reference: bool
        :param calibrate_electric: A boolean flag that specifies whether the electric channels should be calibrated
                                   from units of mV to mv/km. This needs dipole length details. The default is False.
        :type calibrate_electric: bool
        :param calibrate_magnetic: A boolean flag that specifies whether the magnetic channels should be calibrated
                                   to units of mV to nT. The calibration process currently supports Metronix sensor
                                   coils. The default is False.
        :type calibrate_magnetic: bool
        :param calibration_data_electric: A dictionary containing the distances from the center to each electrode.
                                          It should include two sub-dictionaries with the following distances:
                                          - 'ex':
                                              - 'x1': Distance (float) from the center to the North electrode
                                              - 'x2': Distance (float) from the center to the South electrode
                                          - 'ey':
                                              - 'y1': Distance (float) from the center to the East electrode
                                              - 'y2': Distance (float) from the center to the West electrode
        :type calibration_data_electric: dict[str, dict[str, float]]
        :param calibration_data_magnetic: TODO
        :type calibration_data_magnetic: dict
        :param apply_notch_filter: A boolean flag indicating whether a notch filter should be applied. Default is False.
        :type apply_notch_filter: bool
        :param notch_frequency: The frequency to be removed using the notch filter. TODO: Harmonics. Default is None.
        :type notch_frequency: float
        :param process_mt: A boolean flag indicating whether the MT impedance is estimated during processing.
                           Default is True. Set to False to skip MT impedance calculations if only Tipper is
                           being estimated.
        :type process_mt: bool
        :param process_tipper: A boolean flag indicating whether the Tipper is estimated during processing.
                               Default is True. Set to False if Tipper is not part of the survey to avoid unnecessary
                               calculations.

        :return: None
        :rtype: NoneType

        """
        # Set attributes from parameters
        self.sampling_frequency = sampling_frequency
        self.fft_length = fft_length
        self.parzen_window_radius = parzen_window_radius
        self.overlap = overlap
        self.remote_reference = remote_reference
        self.calibration_data_electric = calibration_data_electric
        self.calibration_data_magnetic = calibration_data_magnetic
        self.notch_frequency = notch_frequency
        self.process_mt = process_mt
        self.process_tipper = process_tipper

        # Class related attributes
        self.channels = None
        self.fft_frequencies = None
        self.bandavg_ds = None
        self.avgf = None
        self.dof = None
        self.xfft = None

        # Dividing time series into several time windows of length equals to fft length. Number of windows will
        # depend on the time series overlap.
        self.time_series = _reshape_time_series_with_overlap(time_series=time_series, fft_length=fft_length, overlap=overlap)
        del time_series

        self.ft_list = utils.targetfreq(self.sampling_frequency, self.parzen_window_radius, self.fft_length,
                                       frequencies_per_decade)
        self.get_channels()  # Get list of available ts channel. 'Ex', 'Ey', ....
        if calibrate_electric:
            self.calibrate_electric()
        if apply_notch_filter:
            self.apply_notch()
        self.detrend_time_series()
        self.perform_fft()
        if calibrate_magnetic:
            self.calibrate_mag()
        del self.calibration_data_magnetic
        del self.time_series
        gc.collect()
        self.perform_bandavg()

    def get_channels(self) -> None:
        """
        Get details of channels in the received data.

        :return: None
        :rtype: NoneType

        """

        self.channels = list(self.time_series.keys())

    def calibrate_electric(self) -> None:
        """
        Calibrate electric field data

        :return: None
        :rtype: NoneType

        """
        print('Calibrating electric field channels.')
        if 'ex' in self.channels:
            dipole_ns = abs(self.calibration_data_electric['ex']['x1']) + abs(self.calibration_data_electric['ex']['x2'])
            self.time_series['ex'] = self.time_series['ex'] / (1 * dipole_ns / 1000)
        if 'ey' in self.channels:
            dipole_ew = abs(self.calibration_data_electric['ey']['y1']) + abs(self.calibration_data_electric['ey']['y2'])
            self.time_series['ey'] = self.time_series['ey'] / (1 * dipole_ew / 1000)

    def calibrate_mag(self) -> None:
        """
        Calibrates the magnetic field channels.

        :return: None
        :rtype: NoneType

        """
        print('Calibrating magnetic field channels.')
        desired_elements = ['hx', 'hy', 'hz', 'rx', 'ry']
        # Create a list mag channels out of desired elements if existing in self.channels
        magnetic_channels = [element for element in desired_elements if element in self.channels]
        for channel in magnetic_channels:
            if self.calibration_data_magnetic['instrument'] == 'metronix':
                sensor_type = self.calibration_data_magnetic[channel]['sensor_type']
                sensor_serial_number = str(self.calibration_data_magnetic[channel]['sensor_serial_number'])
                if self.calibration_data_magnetic[channel]['chopper_status'] == 1:
                    chopper_status = "chopper_on"
                elif self.calibration_data_magnetic[channel]['chopper_status'] == 0:
                    chopper_status = "chopper_off"
                else:
                    chopper_status = None
                calibration_data = self.calibration_data_magnetic[channel]['calibration_data'][sensor_serial_number][
                    chopper_status]
                calibration_object = calibration.MetronixCalibration(self.xfft[channel], self.fft_frequencies, sensor_type,
                                                                     chopper_status, calibration_data)
                self.xfft[channel] = calibration_object.calibrated_data
            else:
                raise NotImplementedError(
                    f'Calibration for {self.calibration_data_magnetic["instrument"]} instruments is not implemented as of now.'
                )

    def apply_notch(self) -> None:
        """
        Apply the notch filter

        :return: None
        :rtype: NoneType

        """
        print("Applying notch filter...")
        with ThreadPoolExecutor(max_workers=len(self.time_series)) as executor:
            # Dictionary to hold futures
            futures = {}

            # Submit each task to the executor
            for channel in self.channels:
                futures[channel] = executor.submit(sp.notchfilsos, self.time_series[channel], self.sampling_frequency,
                                                   self.notch_frequency)

            # Retrieve the results when needed
            for channel, future in futures.items():
                self.time_series[channel] = future.result()
        print("Notch filter applied")

    def detrend_time_series(self) -> None:
        """
        Detrend time series

        :return: None
        :rtype: NoneType

        """
        print('Applying detrend to the time series.')
        for channel in self.channels:
            self.time_series[channel] = signal.detrend(self.time_series[channel], axis=0)

    def perform_fft(self) -> None:
        """
        Perform FFT

        :return: None
        :rtype: NoneType

        """
        print('Performing FFT.')
        self.xfft = {}
        for channel in self.channels:
            self.fft_frequencies, self.xfft[channel] = sp.do_fft(self.time_series[channel], self.sampling_frequency,
                                                           self.fft_length)

    def perform_bandavg(self) -> None:
        """
        Perform band averaging.

        :return: None
        :rtype: NoneType

        """
        print('Starting band averaging.')
        self.dof = np.empty(self.ft_list.shape[0], dtype=int)
        self.avgf = np.empty(self.ft_list.shape[0], dtype=int)

        parzen_window = np.empty(
            (self.xfft[next(iter(self.xfft))].shape[0], 1, self.ft_list.shape[0]), dtype=float)

        for i in range(self.ft_list.shape[0]):
            ft = self.ft_list[i]
            parzen_window[:, :, i] = stats.parzen(self.fft_frequencies, ft, self.parzen_window_radius)
            self.dof[i] = (2 * 2 * np.sum(parzen_window[:, :, i] != 0)) - 4
            self.avgf[i] = np.sum(parzen_window[:, :, i] != 0)

        self.bandavg_ds = xr.Dataset(
            coords={
                'time_window': np.arange(self.xfft[next(iter(self.xfft))].shape[1]),
                'frequency': self.ft_list
            }
        )

        if self.process_mt:
            sum_parzen = np.sum(parzen_window, axis=0)

            # Compute the weighted sums
            self.bandavg_ds['ex'] = (
                ('time_window', 'frequency'),
                np.sum(self.xfft['ex'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['ey'] = (
                ('time_window', 'frequency'),
                np.sum(self.xfft['ey'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hx'] = (
                ('time_window', 'frequency'),
                np.sum(self.xfft['hx'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hy'] = (
                ('time_window', 'frequency'),
                np.sum(self.xfft['hy'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            if self.process_tipper:
                self.bandavg_ds['hz'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['hz'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            if self.remote_reference:
                self.bandavg_ds['rx'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['rx'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
                self.bandavg_ds['ry'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['ry'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            # Compute the auto- and cross-spectra
            ex_conj = np.conj(self.xfft['ex'])
            ey_conj = np.conj(self.xfft['ey'])
            if self.process_tipper:
                hz_conj = np.conj(self.xfft['hz'])
            if self.remote_reference:
                hx_conj = np.conj(self.xfft['rx'])
                hy_conj = np.conj(self.xfft['ry'])
            else:
                hx_conj = np.conj(self.xfft['hx'])
                hy_conj = np.conj(self.xfft['hy'])

            self.bandavg_ds['exex'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * ex_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['eyey'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['hxhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hx'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['hyhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hy'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            if self.process_tipper:
                self.bandavg_ds['hzhz'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hz_conj[:, :, np.newaxis] * parzen_window,
                    axis=0) / sum_parzen)
            #
            self.bandavg_ds['exey'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['hxhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hx'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['hyhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hy'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            # Ex output =====
            self.bandavg_ds['exhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['exhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            # Ey output =====
            self.bandavg_ds['eyhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.bandavg_ds['eyhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            if self.process_tipper:
                # Hz output =====
                self.bandavg_ds['hzhx'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                    axis=0) / sum_parzen)
                self.bandavg_ds['hzhy'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                    axis=0) / sum_parzen)

            z_deno = (self.bandavg_ds['hxhx'] * self.bandavg_ds['hyhy']) - (
                    self.bandavg_ds['hxhy'] * self.bandavg_ds['hyhx'])
            #
            zxx_num = (self.bandavg_ds['hyhy'] * self.bandavg_ds['exhx']) - (
                    self.bandavg_ds['hyhx'] * self.bandavg_ds['exhy'])
            zxy_num = (self.bandavg_ds['hxhx'] * self.bandavg_ds['exhy']) - (
                    self.bandavg_ds['hxhy'] * self.bandavg_ds['exhx'])
            self.bandavg_ds['zxx_single'] = zxx_num / z_deno
            self.bandavg_ds['zxy_single'] = zxy_num / z_deno
            #
            zyx_num = (self.bandavg_ds['hyhy'] * self.bandavg_ds['eyhx']) - (
                    self.bandavg_ds['hyhx'] * self.bandavg_ds['eyhy'])
            zyy_num = (self.bandavg_ds['hxhx'] * self.bandavg_ds['eyhy']) - (
                    self.bandavg_ds['hxhy'] * self.bandavg_ds['eyhx'])
            self.bandavg_ds['zyx_single'] = zyx_num / z_deno
            self.bandavg_ds['zyy_single'] = zyy_num / z_deno

            if self.process_tipper:
                t_deno = (self.bandavg_ds['hxhx'] * self.bandavg_ds['hyhy']) - (
                        self.bandavg_ds['hxhy'] * self.bandavg_ds['hyhx'])
                self.bandavg_ds['tzx_single'] = ((self.bandavg_ds['hzhx'] * self.bandavg_ds['hyhy']) - (
                        self.bandavg_ds['hzhy'] * self.bandavg_ds['hyhx'])) / t_deno
                self.bandavg_ds['tzy_single'] = ((self.bandavg_ds['hzhy'] * self.bandavg_ds[
                    'hxhx']) - (self.bandavg_ds['hzhx'] * self.bandavg_ds['hxhy'])) / t_deno

            # Preparing selection arrays

            self.bandavg_ds['ex_selection_coh'] = xr.DataArray(
                np.full(self.bandavg_ds['ex'].shape, True),
                coords=self.bandavg_ds.coords,
                dims=self.bandavg_ds.dims
            )

            self.bandavg_ds['ey_selection_coh'] = xr.DataArray(
                np.full(self.bandavg_ds['ey'].shape, True),
                coords=self.bandavg_ds.coords,
                dims=self.bandavg_ds.dims
            )

            if self.process_tipper:
                self.bandavg_ds['hz_selection_coh'] = xr.DataArray(
                    np.full(self.bandavg_ds['hz'].shape, True),
                    coords=self.bandavg_ds.coords,
                    dims=self.bandavg_ds.dims
                )

            self.bandavg_ds['alpha_e_selection'] = xr.DataArray(
                np.full(self.bandavg_ds['ex'].shape, True),
                coords=self.bandavg_ds.coords,
                dims=self.bandavg_ds.dims
            )

            self.bandavg_ds['alpha_h_selection'] = xr.DataArray(
                np.full(self.bandavg_ds['hx'].shape, True),
                coords=self.bandavg_ds.coords,
                dims=self.bandavg_ds.dims
            )
            print('Band averaging finished.')


def _reshape_time_series_with_overlap(time_series, fft_length, overlap):
    """
    Reshape array with overlap.

    :return: None
    :rtype: NoneType

    """
    for channel in time_series:
        time_series[channel] = utils.reshape_array_with_overlap(window_length=fft_length,
                                                                overlap=overlap,
                                                                data=time_series[channel])

    return time_series
