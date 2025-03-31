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
                 time_series: dict,
                 sampling_frequency: float,
                 fft_length: Optional[int] = 1024,
                 parzen_radius: Optional[float] = 0.25,
                 overlap: Optional[int] = 50,
                 frequencies_per_decade: Optional[int] = 12,
                 remote_reference: Optional[bool] = False,
                 calibrate_magnetic: Optional[bool] = False,
                 calibrate_electric: Optional[bool] = False,
                 calibration_data_electric: Optional[dict] = None,
                 calibration_data_magnetic: Optional[dict] = None,
                 notch_filter_apply: Optional[bool] = False,
                 notch_frequency: Optional[float] = None,
                 process_mt: Optional[bool] = True,
                 process_tipper: Optional[bool] = True) -> None:
        """
        Constructor

        :return: None
        :rtype: NoneType

        """
        # Set attributes from parameters
        self.sampling_frequency = sampling_frequency
        self.fft_length = fft_length
        self.parzen_radius = parzen_radius
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

        self.ft_list = utils.targetfreq(self.sampling_frequency, self.parzen_radius, self.fft_length,
                                       frequencies_per_decade)
        self.get_channels()  # Get list of available ts channel. 'Ex', 'Ey', ....
        if calibrate_electric:
            self.calibrate_electric()
        if notch_filter_apply:
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
            parzen_window[:, :, i] = stats.parzen(self.fft_frequencies, ft, self.parzen_radius)
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
