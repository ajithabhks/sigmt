"""
Class for band averaging. Here, auto and cross spectra are calculated for
each target frequencies and transfer functions are computed for all events.
"""
import gc
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr
from scipy import signal

import sigmt.core.sigproc as sp
import sigmt.core.statistics as stats
import sigmt.utils.utils as utils
from sigmt.utils.metronix.calibration import Calibration


class BandAvg:
    """
    Class to perform band averaging
    """

    def __init__(self, procinfo: dict, bandavg_msg: dict) -> None:
        """
        Constructor
        """
        self.channels = None
        self.fft_freqs = None
        self.bandavg_ds = None
        self.avgf = None
        self.dof = None
        self.data_dict = None
        self.mag_channels = None
        self.xfft = None
        self.procinfo = procinfo
        self.processing_mode = self.procinfo['processing_mode']
        self.fs = procinfo['fs']
        self.fft_length = procinfo['fft_length']
        self.parzen_radius = procinfo['parzen_radius']
        self.notch_status = procinfo['notch']
        self.notch_frequency = procinfo['notch_frequency']
        freq_per_decade = procinfo['frequencies_per_decade']
        #
        self.ts = bandavg_msg['ts']
        self.header = bandavg_msg['header']
        self.cal_data = bandavg_msg['caldata']
        #
        self.ftlist = utils.targetfreq(self.fs, self.parzen_radius, self.fft_length, freq_per_decade)
        self.getchannels()  # Get ts channel info, 'Ex', 'Ey', ....
        self.calibrate_electric()
        if self.notch_status == 'on':
            self.apply_notch()
        self.detrend_ts()
        self.perform_fft()
        self.calibrate_mag()
        del self.cal_data
        del self.ts
        gc.collect()
        self.perform_bandavg()

    def getchannels(self) -> None:
        """
        Get details of channels in the received data.
        """

        self.channels = list(self.ts.keys())

    def calibrate_electric(self) -> None:
        """
        Calibrate electric field data
        """
        print('Calibrating electric field channels.')
        if 'ex' in self.channels:
            dipole_ns = abs(self.header['ex']['x1'][0]) + abs(self.header['ex']['x2'][0])
            self.ts['ex'] = self.ts['ex'] / (1 * dipole_ns / 1000)
        if 'ey' in self.channels:
            dipole_ew = abs(self.header['ey']['y1'][0]) + abs(self.header['ey']['y2'][0])
            self.ts['ey'] = self.ts['ey'] / (1 * dipole_ew / 1000)

    def apply_notch(self) -> None:
        """
        Apply the notch filter
        """
        print("Applying notch filter...")
        with ThreadPoolExecutor(max_workers=len(self.ts)) as executor:
            # Dictionary to hold futures
            futures = {}

            # Submit each task to the executor
            for channel in self.channels:
                futures[channel] = executor.submit(sp.notchfilsos, self.ts[channel], self.fs, self.notch_frequency)

            # Retrieve the results when needed
            for channel, future in futures.items():
                self.ts[channel] = future.result()
        print("Notch filter applied")

    def detrend_ts(self) -> None:
        """
        Detrend time series
        """
        print('Applying detrend to the time series.')
        for channel in self.channels:
            self.ts[channel] = signal.detrend(self.ts[channel], axis=0)

    def perform_fft(self) -> None:
        """
        Perform FFT
        """
        print('Performing FFT.')
        self.xfft = {}
        for channel in self.channels:
            self.fft_freqs, self.xfft[channel] = sp.do_fft(self.ts[channel], self.fs, self.fft_length)

    def calibrate_mag(self) -> None:
        """
        Calibrates the magnetic field channels
        """
        print('Calibrating magnetic field channels.')
        desired_elements = ['hx', 'hy', 'hz', 'rx', 'ry']
        # Create a list mag channels out of desired elements if existing in self.channels
        self.mag_channels = [element for element in desired_elements if element in self.channels]
        for channel in self.mag_channels:
            sensor_type = self.header[channel]['sensor']
            if self.header[channel]['bychopper'][0] == 1:
                stat = "ChoppOn"
            elif self.header[channel]['bychopper'][0] == 0:
                stat = "ChoppOff"
            calibration_data = self.cal_data[str(self.header[channel]['sensor_no'][0])][stat]
            calibration_object = Calibration(self.xfft[channel], self.fft_freqs, sensor_type, stat, calibration_data)
            self.xfft[channel] = calibration_object.calibrated_data

    def perform_bandavg(self) -> None:
        """
        Perform band averaging.
        """
        print('Starting band averaging.')
        self.dof = np.empty(self.ftlist.shape[0], dtype=int)
        self.avgf = np.empty(self.ftlist.shape[0], dtype=int)

        parzen_window = np.empty((self.xfft[next(iter(self.xfft))].shape[0], 1, self.ftlist.shape[0]), dtype=float)

        for i in range(self.ftlist.shape[0]):
            ft = self.ftlist[i]
            parzen_window[:, :, i] = stats.parzen(self.fft_freqs, ft, self.parzen_radius)
            self.dof[i] = (2 * 2 * np.sum(parzen_window[:, :, i] != 0)) - 4
            self.avgf[i] = np.sum(parzen_window[:, :, i] != 0)

        self.bandavg_ds = xr.Dataset(
            coords={
                'time_window': np.arange(self.xfft[next(iter(self.xfft))].shape[1]),
                'frequency': self.ftlist
            }
        )

        if self.processing_mode == "MT + Tipper" or self.processing_mode == "MT Only":
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
            if not self.processing_mode == "MT Only":
                self.bandavg_ds['hz'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['hz'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            if self.procinfo['remotesite'] is not None:
                self.bandavg_ds['rx'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['rx'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
                self.bandavg_ds['ry'] = (
                    ('time_window', 'frequency'),
                    np.sum(self.xfft['ry'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            # Compute the auto- and cross-spectra
            ex_conj = np.conj(self.xfft['ex'])
            ey_conj = np.conj(self.xfft['ey'])
            if not self.processing_mode == "MT Only":
                hz_conj = np.conj(self.xfft['hz'])
            if self.procinfo['remotesite'] is not None:
                hx_conj = np.conj(self.xfft['rx'])
                hy_conj = np.conj(self.xfft['ry'])
            else:
                hx_conj = np.conj(self.xfft['hx'])
                hy_conj = np.conj(self.xfft['hy'])

            self.bandavg_ds['exex'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * ex_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['eyey'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hxhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hx'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hyhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hy'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            if not self.processing_mode == "MT Only":
                self.bandavg_ds['hzhz'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hz_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            #
            self.bandavg_ds['exey'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hxhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hx'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['hyhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['hy'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            # Ex output =====
            self.bandavg_ds['exhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['exhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ex'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            # Ey output =====
            self.bandavg_ds['eyhx'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.bandavg_ds['eyhy'] = (('time_window', 'frequency'), np.sum(
                self.xfft['ey'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            if not self.processing_mode == "MT Only":
                # Hz output =====
                self.bandavg_ds['hzhx'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
                self.bandavg_ds['hzhy'] = (('time_window', 'frequency'), np.sum(
                    self.xfft['hz'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

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

            if not self.processing_mode == "MT Only":
                t_deno = (self.bandavg_ds['hxhx'] * self.bandavg_ds['hyhy']) - (
                        self.bandavg_ds['hxhy'] * self.bandavg_ds['hyhx'])
                self.bandavg_ds['tzx_single'] = ((self.bandavg_ds['hzhx'] * self.bandavg_ds['hyhy']) - (
                        self.bandavg_ds['hzhy'] * self.bandavg_ds['hyhx'])) / t_deno
                self.bandavg_ds['tzy_single'] = ((self.bandavg_ds['hzhy'] * self.bandavg_ds['hxhx']) - (
                        self.bandavg_ds['hzhx'] * self.bandavg_ds['hxhy'])) / t_deno

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

            if not self.processing_mode == "MT Only":
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
