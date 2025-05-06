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
import sigmt.utils.utils as utils
from sigmt.utils.metronix.calibration import MetronixCalibration


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
                 notch_harmonics: Optional[int] = None,
                 process_mt: bool = True,
                 process_tipper: bool = True) -> None:
        """
        Constructor

        TODO: Add docs
            - No separate cross powers for local and remote.
            - Hx and Hy must be given as input.

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
        :param notch_harmonics: Defines the number of harmonics to be filtered using a notch filter.
                                If set to `None` (default), all harmonics up to the highest frequency of interest
                                will be removed. Else, if `notch_frequency` is 50 Hz and `notch_harmonics`
                                is 3, the harmonics at 100 Hz, 150 Hz and 200 Hz will be filtered out.
        :type notch_harmonics: Optional[int]
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
        self.notch_harmonics = notch_harmonics
        self.process_mt = process_mt
        self.process_tipper = process_tipper

        # Class related attributes
        self.channels = None
        self.fft_frequencies = None
        self.band_averaged_dataset = None
        self.avgf = None
        self.dof = None
        self.spectra = None

        # Dividing time series into several time windows of length equals to fft length. Number of windows will
        # depend on the time series overlap.
        self.time_series = _reshape_time_series_with_overlap(time_series=time_series, fft_length=fft_length, overlap=overlap)
        del time_series

        self.ft_list = utils.get_target_frequency_list(sampling_frequency=self.sampling_frequency,
                                                       parzen_window_radius=self.parzen_window_radius,
                                                       fft_length=self.fft_length,
                                                       table_type='default',
                                                       frequencies_per_decade=frequencies_per_decade)

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
        self.perform_band_averaging()

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
                calibration_object = MetronixCalibration(self.spectra[channel], self.fft_frequencies,
                                                         sensor_type,
                                                         chopper_status, calibration_data)
                self.spectra[channel] = calibration_object.calibrated_data
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
        # TODO: Replace this with multiprocessing
        with ThreadPoolExecutor(max_workers=len(self.time_series)) as executor:
            # Dictionary to hold futures
            futures = {}

            # Submit each task to the executor
            for channel in self.channels:
                futures[channel] = executor.submit(sp.notch_filter_sos,
                                                   time_series=self.time_series[channel],
                                                   sampling_frequency=self.sampling_frequency,
                                                   notch_frequency=self.notch_frequency,
                                                   harmonics=self.notch_harmonics)

            # Retrieve the results when needed
            for channel, future in futures.items():
                self.time_series[channel] = future.result()
        print("Notch filter applied.")

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
        self.spectra = {}
        for channel in self.channels:
            self.fft_frequencies, self.spectra[channel] = sp.do_fft(self.time_series[channel], self.sampling_frequency,
                                                                    self.fft_length)

    def perform_band_averaging(self) -> None:
        """
        Perform band averaging.

        :return: None
        :rtype: NoneType

        """
        print('Starting band averaging.')

        # Creating empty arrays
        self.dof = np.empty(self.ft_list.shape[0], dtype=int)  # Degree of freedom
        self.avgf = np.empty(self.ft_list.shape[0], dtype=int)  # Number of frequencies used for averaging

        # Create a 3D array for parzen window for all target frequencies
        parzen_window = np.empty(
            (self.spectra[next(iter(self.spectra))].shape[0], 1, self.ft_list.shape[0]), dtype=float)

        # Populating parzen window arrays based on target frequency and window radius
        for i in range(self.ft_list.shape[0]):
            ft = float(self.ft_list[i])
            parzen_window[:, :, i] = stats.parzen(self.fft_frequencies, ft, self.parzen_window_radius)
            self.dof[i] = (2 * 2 * np.sum(parzen_window[:, :, i] != 0)) - 4
            self.avgf[i] = np.sum(parzen_window[:, :, i] != 0)

        # Create an empty xarray dataset
        self.band_averaged_dataset = xr.Dataset(
            coords={
                'time_window': np.arange(self.spectra[next(iter(self.spectra))].shape[1]),
                'frequency': self.ft_list
            }
        )

        # Precompute the sum of the Parzen window values to optimize performance
        # by avoiding redundant calculations during subsequent operations.
        sum_parzen = np.sum(parzen_window, axis=0)

        # =========== Band averaging ===============================

        # Calculate the cross-power spectra for Hx and Hy, which are common for both MT and Tipper.
        self.band_averaged_dataset['hx'] = (
            ('time_window', 'frequency'),
            np.sum(self.spectra['hx'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
        self.band_averaged_dataset['hy'] = (
            ('time_window', 'frequency'),
            np.sum(self.spectra['hy'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

        if self.remote_reference:
            # If remote referencing is requested, assign the conjugate term to the remote channels.
            self.band_averaged_dataset['rx'] = (
                ('time_window', 'frequency'),
                np.sum(self.spectra['rx'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.band_averaged_dataset['ry'] = (
                ('time_window', 'frequency'),
                np.sum(self.spectra['ry'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            hx_conj = np.conj(self.spectra['rx'])
            hy_conj = np.conj(self.spectra['ry'])
        else:
            # Else, local channels
            hx_conj = np.conj(self.spectra['hx'])
            hy_conj = np.conj(self.spectra['hy'])

        self.band_averaged_dataset['hxhx'] = (('time_window', 'frequency'), np.sum(
            self.spectra['hx'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
            axis=0) / sum_parzen)
        self.band_averaged_dataset['hyhy'] = (('time_window', 'frequency'), np.sum(
            self.spectra['hy'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
            axis=0) / sum_parzen)

        self.band_averaged_dataset['hxhy'] = (('time_window', 'frequency'), np.sum(
            self.spectra['hx'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
            axis=0) / sum_parzen)
        self.band_averaged_dataset['hyhx'] = (('time_window', 'frequency'), np.sum(
            self.spectra['hy'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
            axis=0) / sum_parzen)

        # Common for MT and Tipper
        denominator = ((self.band_averaged_dataset['hxhx'] * self.band_averaged_dataset['hyhy']) -
                       (self.band_averaged_dataset['hxhy'] * self.band_averaged_dataset['hyhx']))

        # TODO: This may not be created in this class
        self.band_averaged_dataset['alpha_h_selection'] = xr.DataArray(
            np.full(self.band_averaged_dataset['hx'].shape, True),
            coords=self.band_averaged_dataset.coords,
            dims=self.band_averaged_dataset.dims
        )

        # Compute if mt channels are requested
        if self.process_mt:
            self.band_averaged_dataset['ex'] = (
                ('time_window', 'frequency'),
                np.sum(self.spectra['ex'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)
            self.band_averaged_dataset['ey'] = (
                ('time_window', 'frequency'),
                np.sum(self.spectra['ey'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            # Compute the auto- and cross-spectra
            ex_conj = np.conj(self.spectra['ex'])
            ey_conj = np.conj(self.spectra['ey'])

            self.band_averaged_dataset['exex'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ex'][:, :, np.newaxis] * ex_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['eyey'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ey'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['exey'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ex'][:, :, np.newaxis] * ey_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            self.band_averaged_dataset['exhx'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ex'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['exhy'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ex'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            self.band_averaged_dataset['eyhx'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ey'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['eyhy'] = (('time_window', 'frequency'), np.sum(
                self.spectra['ey'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            # Computing the MT impedance for all time windows
            zxx_num = ((self.band_averaged_dataset['hyhy'] * self.band_averaged_dataset['exhx']) -
                       (self.band_averaged_dataset['hyhx'] * self.band_averaged_dataset['exhy']))
            zxy_num = ((self.band_averaged_dataset['hxhx'] * self.band_averaged_dataset['exhy']) -
                       (self.band_averaged_dataset['hxhy'] * self.band_averaged_dataset['exhx']))
            self.band_averaged_dataset['zxx'] = zxx_num / denominator
            self.band_averaged_dataset['zxy'] = zxy_num / denominator
            #
            zyx_num = ((self.band_averaged_dataset['hyhy'] * self.band_averaged_dataset['eyhx']) -
                       (self.band_averaged_dataset['hyhx'] * self.band_averaged_dataset['eyhy']))
            zyy_num = ((self.band_averaged_dataset['hxhx'] * self.band_averaged_dataset['eyhy']) -
                       (self.band_averaged_dataset['hxhy'] * self.band_averaged_dataset['eyhx']))
            self.band_averaged_dataset['zyx'] = zyx_num / denominator
            self.band_averaged_dataset['zyy'] = zyy_num / denominator

            # Preparing selection arrays, this may be used for data selection later.
            # TODO: This may not be created in this class
            self.band_averaged_dataset['ex_selection_coh'] = xr.DataArray(
                np.full(self.band_averaged_dataset['ex'].shape, True),
                coords=self.band_averaged_dataset.coords,
                dims=self.band_averaged_dataset.dims
            )
            self.band_averaged_dataset['ey_selection_coh'] = xr.DataArray(
                np.full(self.band_averaged_dataset['ey'].shape, True),
                coords=self.band_averaged_dataset.coords,
                dims=self.band_averaged_dataset.dims
            )

            # TODO: This may not be created in this class
            self.band_averaged_dataset['alpha_e_selection'] = xr.DataArray(
                np.full(self.band_averaged_dataset['ex'].shape, True),
                coords=self.band_averaged_dataset.coords,
                dims=self.band_averaged_dataset.dims
            )

        # Compute if tipper channels are requested
        if self.process_tipper:
            self.band_averaged_dataset['hz'] = (
                ('time_window', 'frequency'),
                np.sum(self.spectra['hz'][:, :, np.newaxis] * parzen_window, axis=0) / sum_parzen)

            hz_conj = np.conj(self.spectra['hz'])

            self.band_averaged_dataset['hzhz'] = (('time_window', 'frequency'), np.sum(
                self.spectra['hz'][:, :, np.newaxis] * hz_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['hzhx'] = (('time_window', 'frequency'), np.sum(
                self.spectra['hz'][:, :, np.newaxis] * hx_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)
            self.band_averaged_dataset['hzhy'] = (('time_window', 'frequency'), np.sum(
                self.spectra['hz'][:, :, np.newaxis] * hy_conj[:, :, np.newaxis] * parzen_window,
                axis=0) / sum_parzen)

            # Computing the Tipper transfer function for all time windows
            self.band_averaged_dataset['tzx'] = ((self.band_averaged_dataset['hzhx'] * self.band_averaged_dataset[
                'hyhy']) - (self.band_averaged_dataset['hzhy'] * self.band_averaged_dataset['hyhx'])) / denominator
            self.band_averaged_dataset['tzy'] = ((self.band_averaged_dataset['hzhy'] * self.band_averaged_dataset[
                'hxhx']) - (self.band_averaged_dataset['hzhx'] * self.band_averaged_dataset['hxhy'])) / denominator

            # TODO: This may not be created in this class
            self.band_averaged_dataset['hz_selection_coh'] = xr.DataArray(
                np.full(self.band_averaged_dataset['hz'].shape, True),
                coords=self.band_averaged_dataset.coords,
                dims=self.band_averaged_dataset.dims
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
