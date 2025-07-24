"""
Class for the robust estimation of MT impedance and tipper
"""
import time

import numpy as np
import xarray as xr
from sklearn.covariance import MinCovDet

import sigmt.core.statistics as stats


class RobustEstimation:
    """
    Class for robust estimation
    """

    def __init__(self, procinfo: dict, dataset: xr.Dataset) -> None:
        """
        Constructor

        :param procinfo: Dictionary containing information for processing such as target
                        frequency, local and remote station details, FFT length, parzen radius,
                        location, etc.
                        :TODO make it instrument independent
        :type procinfo: dict
        :param dataset: xarray dataset containing required numpy arrays for the regression.
                        Cross powers, Selection booleans, TF estimates for all time windows.

        :return: None
        :rtype: NoneType

        """
        self.ft = None
        self.output_channels = None
        self.channel = None
        self.filtered_dataset = None
        self.estimates = None
        self.z2_initial_jackknife = None
        self.z1_initial_jackknife = None
        self.z2_mean_md = None
        self.z1_mean_md = None
        self.output = None
        self.z1_num_keys = None
        self.z2_num_keys = None
        self.z_deno_keys = None
        self.residuals = None
        self.huber_weights = None
        self.z1_robust_huber = None
        self.z2_robust_huber = None
        self.z1_var = None
        self.z2_var = None
        self.predicted_coherency = None
        self.procinfo = procinfo
        self.processing_mode = self.procinfo['processing_mode']
        self.get_output_channels()
        self.dataset = dataset
        self.target_frequencies = dataset['frequency'].values
        self.get_estimate_and_variance()

    def get_output_channels(self) -> None:
        """
        Creates output channels based on processing mode.
        SigMT can support MT+Tipper, MT Only and Tipper Only mode.

        :return: None
        :rtype: NoneType

        """
        if self.processing_mode == "MT + Tipper":
            self.output_channels = ['ex', 'ey', 'hz']
        elif self.processing_mode == "MT Only":
            self.output_channels = ['ex', 'ey']
        elif self.processing_mode == "Tipper Only":
            self.output_channels = ['hz']

    def get_estimate_and_variance(self) -> None:
        """
        Computes the robust estimates and variances.
        This method loops over available output channels and target frequencies.
        It rejects noisy windows and pass data to robust estimation and variance
        calculation methods and saves data in an xarray.

        :return: None
        :rtype: NoneType

        """
        estimate = {}
        estimation_time = time.time()
        for self.channel in self.output_channels:
            z1 = []
            z2 = []
            z1_var = []
            z2_var = []
            coh = []
            for self.ft in self.target_frequencies:
                single_time = time.time()
                # Filtering based on target frequency
                self.filtered_dataset = self.dataset.sel(frequency=self.ft)
                # Filtering based on a selection array for ex/ey/hz
                self.get_selection_array()
                selected_windows = \
                    np.where(self.filtered_dataset[f'selection_array_{self.channel}'].values.astype(
                        bool))[0]
                print(
                    f'Channel: {self.channel}, ft: {self.ft} Hz. Applying data selection conditions.')
                self.filtered_dataset = self.filtered_dataset.isel(time_window=selected_windows)
                # Filtering based on mahalanobis distance
                self.get_mahalanobis_distance()
                mahalanobis_selection = self.filtered_dataset['maha_dist'] < self.procinfo[
                    'md_thresh']
                selected_windows = np.where(mahalanobis_selection.astype(bool))[0]
                if len(selected_windows) > 10:
                    print(
                        f'Channel: {self.channel}, ft: {self.ft} Hz. Applying Mahalanobis distance condition.')
                    self.filtered_dataset = self.filtered_dataset.isel(time_window=selected_windows)
                else:
                    print('Skipping Mahalanobis distance condition due to insufficient windows.')
                print(
                    f'Channel: {self.channel}, ft: {self.ft} Hz. Getting Jackknife initial guess.')
                self.get_jackknife_initial_guess()
                print(f'Channel: {self.channel}, ft: {self.ft} Hz. Performing robust estimation.')
                self.perform_robust_estimation()
                print(f'Channel: {self.channel}, ft: {self.ft} Hz. Computing variance.')
                self.get_variance()
                z1.append(self.z1_robust_huber)
                z2.append(self.z2_robust_huber)
                z1_var.append(self.z1_var)
                z2_var.append(self.z2_var)
                coh.append(self.predicted_coherency)
                print(f'Time taken: {time.time() - single_time} seconds.')
            #
            if self.channel == 'ex':
                z1_key = "zxx"
                z2_key = "zxy"
                z1_var_key = "zxx_var"
                z2_var_key = "zxy_var"
                coh_key = 'coh_ex'
            if self.channel == 'ey':
                z1_key = "zyx"
                z2_key = "zyy"
                z1_var_key = "zyx_var"
                z2_var_key = "zyy_var"
                coh_key = 'coh_ey'
            if self.channel == 'hz':
                z1_key = "tzx"
                z2_key = "tzy"
                z1_var_key = "tzx_var"
                z2_var_key = "tzy_var"
                coh_key = 'coh_hz'
            estimate[z1_key] = xr.DataArray(np.asarray(z1),
                                            coords={'frequency': self.target_frequencies},
                                            dims='frequency'
                                            )
            estimate[z2_key] = xr.DataArray(np.asarray(z2),
                                            coords={'frequency': self.target_frequencies},
                                            dims='frequency'
                                            )
            estimate[z1_var_key] = xr.DataArray(np.asarray(z1_var),
                                                coords={'frequency': self.target_frequencies},
                                                dims='frequency'
                                                )
            estimate[z2_var_key] = xr.DataArray(np.asarray(z2_var),
                                                coords={'frequency': self.target_frequencies},
                                                dims='frequency'
                                                )
            estimate[coh_key] = xr.DataArray(np.asarray(coh),
                                             coords={'frequency': self.target_frequencies},
                                             dims='frequency'
                                             )
            self.estimates = xr.Dataset(estimate)
        print(f'Total time taken for robust estimation: {time.time() - estimation_time} seconds.')

    def get_selection_array(self) -> None:
        """
        Prepare an array of boolean suggesting selected time windows.

        :return: None
        :rtype: NoneType

        """
        if not self.processing_mode == "Tipper Only":
            self.filtered_dataset['selection_array_ex'] = (
                    self.filtered_dataset['ex_selection_coh'] * self.filtered_dataset['alpha_e_selection'] *
                    self.filtered_dataset['alpha_h_selection'])
            self.filtered_dataset['selection_array_ey'] = (
                    self.filtered_dataset['ey_selection_coh'] * self.filtered_dataset['alpha_e_selection'] *
                    self.filtered_dataset['alpha_h_selection'])
        if not self.processing_mode == "MT Only":
            self.filtered_dataset['selection_array_hz'] = (
                    self.filtered_dataset['hz_selection_coh'] * self.filtered_dataset['alpha_h_selection'])

        # Avoid crashing due to insufficient data
        if not self.processing_mode == "Tipper Only":
            if np.sum(self.filtered_dataset['selection_array_ex']) < 10:
                print('There is not enough data to continue the regression.')
                print(f'Trying to remove polarization direction related data rejections for ft:{self.ft} Hz (ex) if applied.')
                self.filtered_dataset['selection_array_ex'] = self.filtered_dataset['ex_selection_coh']
                if np.sum(self.filtered_dataset['selection_array_ex']) <= 10:
                    print('There is still not enough data to continue the regression.')
                    print(f'Trying to remove coherency threshold related data rejection for ft:{self.ft} Hz (ex) if applied.')
                    self.filtered_dataset['selection_array_ex'] = xr.ones_like(self.filtered_dataset['selection_array_ex'],
                                                                               dtype=bool)
        if not self.processing_mode == "Tipper Only":
            if np.sum(self.filtered_dataset['selection_array_ey']) < 10:
                print('There is not enough data to continue the regression.')
                print(f'Trying to remove polarization direction related data rejections for ft:{self.ft} Hz (ey) if applied.')
                self.filtered_dataset['selection_array_ey'] = self.filtered_dataset['ey_selection_coh']
                if np.sum(self.filtered_dataset['selection_array_ey']) <= 10:
                    print('There is still not enough data to continue the regression.')
                    print(f'Trying to remove coherency threshold related data rejection for ft:{self.ft} Hz (ey) if applied.')
                    self.filtered_dataset['selection_array_ey'] = xr.ones_like(self.filtered_dataset['selection_array_ey'],
                                                                               dtype=bool)

        if not self.processing_mode == "MT Only":
            if np.sum(self.filtered_dataset['selection_array_hz']) < 10:
                print('There is not enough data to continue the regression.')
                print(
                    f'Trying to remove polarization direction related data rejections for ft:{self.ft} Hz (hz) if applied.')
                self.filtered_dataset['selection_array_hz'] = self.filtered_dataset['hz_selection_coh']
                if np.sum(self.filtered_dataset['selection_array_hz']) <= 10:
                    print('There is still not enough data to continue the regression.')
                    print(
                        f'Trying to remove coherency threshold related data rejection for ft:{self.ft} Hz (hz) if applied.')
                    self.filtered_dataset['selection_array_hz'] = xr.ones_like(self.filtered_dataset['selection_array_hz'],
                                                                               dtype=bool)

    def get_mahalanobis_distance(self) -> None:
        """
        Method to calculate the mahalanobis distance

        :return: None
        :rtype: NoneType

        """
        if self.channel == 'ex':
            z1 = self.filtered_dataset['zxx'].values
            z2 = self.filtered_dataset['zxy'].values
        if self.channel == 'ey':
            z1 = self.filtered_dataset['zyx'].values
            z2 = self.filtered_dataset['zyy'].values
        if self.channel == 'hz':
            z1 = self.filtered_dataset['tzx'].values
            z2 = self.filtered_dataset['tzy'].values
        data = np.transpose(np.vstack((z1.real, z1.imag, z2.real, z2.imag)))
        robust_cov = MinCovDet(random_state=0).fit(data)
        self.filtered_dataset['maha_dist'] = xr.DataArray(
            np.sqrt(robust_cov.mahalanobis(data)),
            coords=self.filtered_dataset.coords,
            dims=self.filtered_dataset.dims
        )
        self.z1_mean_md = complex(robust_cov.location_[0], robust_cov.location_[1])
        self.z2_mean_md = complex(robust_cov.location_[2], robust_cov.location_[3])

    def get_jackknife_initial_guess(self) -> None:
        """
        Prepares jackknife initial guess from transfer function values from selected time windows.

        :return: None
        :rtype: NoneType

        """
        if self.channel == 'ex':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['zxx'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['zxy'])
        if self.channel == 'ey':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['zyx'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['zyy'])
        if self.channel == 'hz':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['tzx'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['tzy'])

    def perform_robust_estimation(self) -> None:
        """
        Performs robust estimation.

        The robust estimation algorithm used here is adapted from the following sources:

        1.  Ritter, O., Junge, A., Dawes, G.J., 1998. New equipment and processing for
            magnetotelluric remote reference observations. Geophys. J. Int. 132 (3), 535–548.
            https://doi.org/10.1046/j.1365-246X.1998.00440.x.

        2.  Manoj, C., 2003. Magnetotelluric Data Analysis Using Advances in Signal Processing
            Techniques. Ph.D. Dissertation. CSIR – National Geophysical Research Institute and
            Osmania University, India.

        3.  Chave, A., Jones, A. (Eds.), 2012. The Magnetotelluric Method: Theory and Practice.
            Cambridge University Press, Cambridge. Pages: 188-189

        :return: None
        :rtype: NoneType

        """
        self.z1_robust_huber = None
        self.z2_robust_huber = None
        prev_sum_squared = None

        self.get_keys()
        hx = self.filtered_dataset['hx']
        hy = self.filtered_dataset['hy']
        n_time_windows = self.output.shape[0]

        z1 = self.z1_initial_jackknife
        z2 = self.z2_initial_jackknife

        # ------------ Iteration starts here ------------------
        for iteration in range(20):
            # Calculating residuals
            self.residuals = abs(self.output - (z1 * hx) - (z2 * hy))

            if iteration == 0:
                # Iteration = 0 produces the preliminary estimates of the transfer function based on MAD.
                # Robust estimation of transfer function (A1.3 - Ritter (1998)) begins with iteration = 1.
                #
                # Initial guess of variance based on MAD scale estimate
                dm = 1.483 * np.median(abs(self.residuals - np.median(self.residuals)))
                # Upper limit to the MAD scale estimate
                scale_factor = 1.5
                km = scale_factor * dm
                # Allowing up to 5% of the data by adjusting scale_factor = 1.5.
                # Because, in some cases, high residual values prevent Huber weight
                # conditions from being met, causing the runtime error when Lc becomes zero.
                while int(np.sum(self.residuals <= km)) < int(np.ceil(len(self.residuals) * 0.05)):
                    i = i + 1
                    scale_factor = scale_factor + 0.1
                    km = scale_factor * dm
                # Get huber weights based on km
                self.get_huber_weights(km)
            else:
                # Number of windows in which weight = 1
                lc = np.sum((self.huber_weights == 1) * 1)
                # Weighted sum of squared residuals
                sum_squared = (n_time_windows / (lc ** 2)) * (
                    np.sum(self.huber_weights * (self.residuals ** 2)))
                if prev_sum_squared is not None:
                    change = abs(prev_sum_squared - sum_squared) / prev_sum_squared
                    if change < 0.01:  # 1%
                        # Stop iteration when there is no change in sum squared more than 1%
                        # This logic is taken from
                        # Chave, A., Jones, A. (Eds.), 2012. The Magnetotelluric Method:
                        # Theory and Practice. Cambridge University Press, Cambridge.
                        # Page: 189
                        break
                prev_sum_squared = sum_squared
                # New variance of the robust solution
                dh = np.sqrt(sum_squared)
                # Upper limit
                kh = float(1.5 * dh)
                # Get huber weights based on kh
                self.get_huber_weights(kh)

            # Applying weights to band averaged cross-spectra
            element_dict = {}
            element_avg_dict = {}
            for i in range(4):
                element_dict[f'z1_num_{i}'] = self.filtered_dataset[
                                                  self.z1_num_keys[i]].values * self.huber_weights
                element_dict[f'z2_num_{i}'] = self.filtered_dataset[
                                                  self.z2_num_keys[i]].values * self.huber_weights
                element_dict[f'z_deno_{i}'] = self.filtered_dataset[
                                                  self.z_deno_keys[i]].values * self.huber_weights
            # Weighted mean of cross-spectra
            for i in range(4):
                element_avg_dict[f'z1_num_avg_{i}'] = np.sum(element_dict[f'z1_num_{i}']) / np.sum(
                    self.huber_weights)
                element_avg_dict[f'z2_num_avg_{i}'] = np.sum(element_dict[f'z2_num_{i}']) / np.sum(
                    self.huber_weights)
                element_avg_dict[f'z_deno_avg_{i}'] = np.sum(element_dict[f'z_deno_{i}']) / np.sum(
                    self.huber_weights)

            #
            z_deno = ((element_avg_dict['z_deno_avg_0'] * element_avg_dict['z_deno_avg_1']) -
                      (element_avg_dict['z_deno_avg_2'] * element_avg_dict['z_deno_avg_3']))
            z1 = (((element_avg_dict['z1_num_avg_0'] * element_avg_dict['z1_num_avg_1']) -
                   (element_avg_dict['z1_num_avg_2'] * element_avg_dict['z1_num_avg_3'])) /
                  z_deno)
            z2 = (((element_avg_dict['z2_num_avg_0'] * element_avg_dict['z2_num_avg_1']) -
                   (element_avg_dict['z2_num_avg_2'] * element_avg_dict['z2_num_avg_3'])) /
                  z_deno)

        if self.channel == 'ex' or self.channel == 'ey':
            # Metronix specific correction
            self.z1_robust_huber = z1 * -1
            self.z2_robust_huber = z2 * -1
        else:
            self.z1_robust_huber = z1
            self.z2_robust_huber = z2

    def get_keys(self) -> None:
        """
        Prepares keys based on the current channel of processing.

        :return: None
        :rtype: NoneType

        """
        self.z_deno_keys = ["hxhx", "hyhy", "hxhy", "hyhx"]
        if self.channel == 'ex':
            self.output = self.filtered_dataset['ex']
            self.z1_num_keys = ["hyhy", "exhx", "hyhx", "exhy"]  # Zxx
            self.z2_num_keys = ["hxhx", "exhy", "hxhy", "exhx"]  # Zxy
        if self.channel == 'ey':
            self.output = self.filtered_dataset['ey']
            self.z1_num_keys = ["hyhy", "eyhx", "hyhx", "eyhy"]  # Zyx
            self.z2_num_keys = ["hxhx", "eyhy", "hxhy", "eyhx"]  # Zyy
        if self.channel == 'hz':
            self.output = self.filtered_dataset['hz']
            self.z1_num_keys = ["hyhy", "hzhx", "hyhx", "hzhy"]  # Zyx
            self.z2_num_keys = ["hxhx", "hzhy", "hxhy", "hzhx"]  # Zyy

    def get_huber_weights(self, upper_limit: float) -> None:
        """
        Prepares huber weights based on the upper limit.
        """
        huber_matrix1 = (self.residuals <= upper_limit) * 1
        huber_matrix2 = (self.residuals > upper_limit) * 1
        huber_matrix2 = huber_matrix2 * (upper_limit / self.residuals)
        self.huber_weights = (huber_matrix1 + huber_matrix2).values

    def get_variance(self) -> None:
        """
        Parametric estimation of variances. Based on the predicted coherency and degree of freedom.
        Dr. Manoj C. Nair helped with this computation.

        :return: None
        :rtype: NoneType

        """
        if self.channel == 'ex':
            outout = np.sum(self.filtered_dataset['exex'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhx = np.sum(self.filtered_dataset['exhx'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhy = np.sum(self.filtered_dataset['exhy'].values * self.huber_weights) / np.sum(
                self.huber_weights)
        if self.channel == 'ey':
            outout = np.sum(self.filtered_dataset['eyey'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhx = np.sum(self.filtered_dataset['eyhx'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhy = np.sum(self.filtered_dataset['eyhy'].values * self.huber_weights) / np.sum(
                self.huber_weights)
        if self.channel == 'hz':
            outout = np.sum(self.filtered_dataset['hzhz'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhx = np.sum(self.filtered_dataset['hzhx'].values * self.huber_weights) / np.sum(
                self.huber_weights)
            outhy = np.sum(self.filtered_dataset['hzhy'].values * self.huber_weights) / np.sum(
                self.huber_weights)

        hxhx = np.sum(self.filtered_dataset['hxhx'].values * self.huber_weights) / np.sum(
            self.huber_weights)
        hxhy = np.sum(self.filtered_dataset['hxhy'].values * self.huber_weights) / np.sum(
            self.huber_weights)
        hyhx = np.sum(self.filtered_dataset['hyhx'].values * self.huber_weights) / np.sum(
            self.huber_weights)
        hyhy = np.sum(self.filtered_dataset['hyhy'].values * self.huber_weights) / np.sum(
            self.huber_weights)

        zpz = self.z1_robust_huber * np.conj(outhx) + self.z2_robust_huber * np.conj(outhy)
        zpx = self.z1_robust_huber * hxhx + self.z2_robust_huber * hyhx
        zpy = self.z1_robust_huber * hxhy + self.z2_robust_huber * hyhy
        zpzp = self.z1_robust_huber * np.conj(zpx) + self.z2_robust_huber * np.conj(zpy)
        zpzp = zpzp * outout

        coh = complex(np.where(np.abs(zpzp) > 0, zpz / np.sqrt(zpzp), 1 + 1j))
        coh = np.abs(coh)

        coh = float(np.where(coh > 1.0, 1 / coh, coh))

        coh = coh ** 2

        dof = int(self.filtered_dataset['dof'].values)

        fis = stats.fisher_critical(dof)

        d = (4 / dof) * fis * (1.0 - coh) * outout

        if abs(hxhx) * abs(hyhy) > 0:
            r2xy = (hxhy * hyhx) / (hxhx * hyhy)
        else:
            r2xy = 100

        if r2xy > 0.999:
            r2xy = 0.999

        if abs(hxhx * (1 - r2xy)) > 0:
            da = d / (hxhx * (1 - r2xy))
        else:
            da = 100

        if abs(hyhy * (1 - r2xy)) > 0:
            db = d / (hyhy * (1 - r2xy))
        else:
            db = 100

        self.z1_var = np.sqrt(abs(np.real(da)))
        self.z2_var = np.sqrt(abs(np.real(db)))
        #
        self.predicted_coherency = np.sqrt(coh)
