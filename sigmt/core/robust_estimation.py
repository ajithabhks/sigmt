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
        self.hx = None
        self.hy = None
        self.z1_num_keys = None
        self.z2_num_keys = None
        self.z_deno_keys = None
        self.residuals = None
        self.kmx = None
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
        Creates output channels based on processing mode
        """
        if self.processing_mode == "MT + Tipper":
            self.output_channels = ['ex', 'ey', 'hz']
        elif self.processing_mode == "MT Only":
            self.output_channels = ['ex', 'ey']
        elif self.processing_mode == "Tipper Only":
            self.output_channels = ['hz']

    def get_estimate_and_variance(self) -> None:
        """
        Computes the robust estimates and variances
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
                    np.where(self.filtered_dataset[f'selection_array_{self.channel}'].values.astype(bool))[0]
                print(f'Channel: {self.channel}, ft: {self.ft} Hz. Applying data selection conditions.')
                self.filtered_dataset = self.filtered_dataset.isel(time_window=selected_windows)
                # Filtering based on mahalanobis distance
                self.get_mahalanobis_distance()
                mahalanobis_selection = self.filtered_dataset['maha_dist'] < self.procinfo['md_thresh']
                selected_windows = np.where(mahalanobis_selection.astype(bool))[0]
                print(f'Channel: {self.channel}, ft: {self.ft} Hz. Applying Mahalanobis distance condition.')
                self.filtered_dataset = self.filtered_dataset.isel(time_window=selected_windows)
                print(f'Channel: {self.channel}, ft: {self.ft} Hz. Getting Jackknife initial guess.')
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
        """
        self.filtered_dataset['selection_array_ex'] = (
                self.filtered_dataset['ex_selection_coh'] * self.filtered_dataset['alpha_e_selection'] *
                self.filtered_dataset['alpha_h_selection'])
        self.filtered_dataset['selection_array_ey'] = (
                self.filtered_dataset['ey_selection_coh'] * self.filtered_dataset['alpha_e_selection'] *
                self.filtered_dataset['alpha_h_selection'])
        if not self.processing_mode == "MT Only":
            self.filtered_dataset['selection_array_hz'] = (
                    self.filtered_dataset['hz_selection_coh'] * self.filtered_dataset['alpha_e_selection'] *
                    self.filtered_dataset['alpha_h_selection'])
        # Avoiding GUI from crashing due to insufficient data
        if np.sum(self.filtered_dataset['selection_array_ex']) < 10:
            print(
                f'Skipping polarization direction selection for ft:{self.ft} Hz (ex) '
                f'due to insufficient data for robust regression. Try with different min and max values.')
            self.filtered_dataset['selection_array_ex'] = self.filtered_dataset['ex_selection_coh']
        #
        if np.sum(self.filtered_dataset['selection_array_ey']) < 10:
            print(
                f'Skipping polarization direction selection for ft:{self.ft} Hz (ey) '
                f'due to insufficient data for robust regression. Try with different min and max values.')
            self.filtered_dataset['selection_array_ey'] = self.filtered_dataset['ey_selection_coh']
        #
        if not self.processing_mode == "MT Only":
            if np.sum(self.filtered_dataset['selection_array_hz']) < 10:
                print(
                    f'Skipping polarization direction selection for ft:{self.ft} Hz (hz) '
                    f'due to insufficient data for robust regression. Try with different min and max values.')
                self.filtered_dataset['selection_array_hz'] = self.filtered_dataset['hz_selection_coh']

    def get_mahalanobis_distance(self) -> None:
        """
        Calculating mahalanobis distance
        """
        if self.channel == 'ex':
            z1 = self.filtered_dataset['zxx_single'].values
            z2 = self.filtered_dataset['zxy_single'].values
        if self.channel == 'ey':
            z1 = self.filtered_dataset['zyx_single'].values
            z2 = self.filtered_dataset['zyy_single'].values
        if self.channel == 'hz':
            z1 = self.filtered_dataset['tzx_single'].values
            z2 = self.filtered_dataset['tzy_single'].values
        data = np.transpose(np.vstack((z1.real, z1.imag, z2.real, z2.imag)))
        robust_cov = MinCovDet().fit(data)
        self.filtered_dataset['maha_dist'] = xr.DataArray(
            np.sqrt(robust_cov.mahalanobis(data)),
            coords=self.filtered_dataset.coords,
            dims=self.filtered_dataset.dims
        )
        self.z1_mean_md = complex(robust_cov.location_[0], robust_cov.location_[1])
        self.z2_mean_md = complex(robust_cov.location_[2], robust_cov.location_[3])

    def get_jackknife_initial_guess(self) -> None:
        """
        Prepares jackknife intial guess
        """
        if self.channel == 'ex':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['zxx_single'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['zxy_single'])
        if self.channel == 'ey':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['zyx_single'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['zyy_single'])
        if self.channel == 'hz':
            self.z1_initial_jackknife = stats.jackknife(self.filtered_dataset['tzx_single'])
            self.z2_initial_jackknife = stats.jackknife(self.filtered_dataset['tzy_single'])

    def perform_robust_estimation(self) -> None:
        """
        Performs robust estimation
        """
        self.get_keys()
        self.hx = self.filtered_dataset['hxhx']
        self.hy = self.filtered_dataset['hyhy']
        n_time_windows = self.output.shape[0]
        # Calculating residuals
        self.residuals = abs(abs(self.output) - (abs(self.z1_initial_jackknife) *
                                                 abs(self.hx)) - (abs(self.z2_initial_jackknife) * abs(self.hy)))
        # Initial variance based on MAD scale estimate
        dmx = 1.483 * np.median(abs(self.residuals - np.median(self.residuals)))
        # Upper limit to the MAD scale estimate
        self.kmx = 1.5 * dmx
        self.get_huber_weights()
        #
        element_dict = {}
        element_avg_dict = {}
        # Applying weights to band averaged cross-spectra
        for i in range(4):
            element_dict[f'z1_num_{i}'] = self.filtered_dataset[self.z1_num_keys[i]].values * self.huber_weights
            element_dict[f'z2_num_{i}'] = self.filtered_dataset[self.z2_num_keys[i]].values * self.huber_weights
            element_dict[f'z_deno_{i}'] = self.filtered_dataset[self.z_deno_keys[i]].values * self.huber_weights
        # Weighted mean of cross-spectra
        for i in range(4):
            element_avg_dict[f'z1_num_avg_{i}'] = np.sum(element_dict[f'z1_num_{i}']) / np.sum(self.huber_weights)
            element_avg_dict[f'z2_num_avg_{i}'] = np.sum(element_dict[f'z2_num_{i}']) / np.sum(self.huber_weights)
            element_avg_dict[f'z_deno_avg_{i}'] = np.sum(element_dict[f'z_deno_{i}']) / np.sum(self.huber_weights)

        z_deno = ((element_avg_dict['z_deno_avg_0'] * element_avg_dict['z_deno_avg_1']) -
                  (element_avg_dict['z_deno_avg_2'] * element_avg_dict['z_deno_avg_3']))
        self.z1_robust_huber = (((element_avg_dict['z1_num_avg_0'] * element_avg_dict['z1_num_avg_1']) -
                                 (element_avg_dict['z1_num_avg_2'] * element_avg_dict['z1_num_avg_3'])) /
                                z_deno)
        self.z2_robust_huber = (((element_avg_dict['z2_num_avg_0'] * element_avg_dict['z2_num_avg_1']) -
                                 (element_avg_dict['z2_num_avg_2'] * element_avg_dict['z2_num_avg_3'])) /
                                z_deno)
        # ------------ Iteration starts here ------------------
        for iteration in range(20):
            lc = np.sum((self.huber_weights == 1) * 1)
            if lc == 0:
                lc = 1
            self.output = self.output * self.huber_weights
            self.hx = self.hx * self.huber_weights
            self.hy = self.hy * self.huber_weights
            self.residuals = abs(abs(self.output) -
                                 (abs(self.z1_robust_huber) * abs(self.hx)) -
                                 (abs(self.z2_robust_huber) * abs(self.hy)))
            # New variance of the robust solution
            dhx = np.sqrt((n_time_windows / (lc ** 2)) * (np.sum(self.huber_weights * (self.residuals ** 2))))
            # Upper limit
            self.kmx = 1.5 * dhx
            self.get_huber_weights()
            #
            # Applying weights to band averaged cross-spectra
            for i in range(4):
                element_dict[f'z1_num_{i}'] = element_dict[f'z1_num_{i}'] * self.huber_weights
                element_dict[f'z2_num_{i}'] = element_dict[f'z2_num_{i}'] * self.huber_weights
                element_dict[f'z_deno_{i}'] = element_dict[f'z_deno_{i}'] * self.huber_weights
            # Weighted mean of cross-spectra
            for i in range(4):
                element_avg_dict[f'z1_num_avg_{i}'] = np.sum(element_dict[f'z1_num_{i}']) / np.sum(self.huber_weights)
                element_avg_dict[f'z2_num_avg_{i}'] = np.sum(element_dict[f'z2_num_{i}']) / np.sum(self.huber_weights)
                element_avg_dict[f'z_deno_avg_{i}'] = np.sum(element_dict[f'z_deno_{i}']) / np.sum(self.huber_weights)

            #
            z_deno = ((element_avg_dict['z_deno_avg_0'] * element_avg_dict['z_deno_avg_1']) -
                      (element_avg_dict['z_deno_avg_2'] * element_avg_dict['z_deno_avg_3']))
            self.z1_robust_huber = (((element_avg_dict['z1_num_avg_0'] * element_avg_dict['z1_num_avg_1']) -
                                     (element_avg_dict['z1_num_avg_2'] * element_avg_dict['z1_num_avg_3'])) /
                                    z_deno)
            self.z2_robust_huber = (((element_avg_dict['z2_num_avg_0'] * element_avg_dict['z2_num_avg_1']) -
                                     (element_avg_dict['z2_num_avg_2'] * element_avg_dict['z2_num_avg_3'])) /
                                    z_deno)
            if self.channel != 'hz':
                self.z1_robust_huber = self.z1_robust_huber * -1
                self.z2_robust_huber = self.z2_robust_huber * -1

    def get_keys(self) -> None:
        """
        Prepares keys
        """
        self.z_deno_keys = ["hxhx", "hyhy", "hxhy", "hyhx"]
        if self.channel == 'ex':
            self.output = self.filtered_dataset['exex']
            self.z1_num_keys = ["hyhy", "exhx", "hyhx", "exhy"]  # Zxx
            self.z2_num_keys = ["hxhx", "exhy", "hxhy", "exhx"]  # Zxy
        if self.channel == 'ey':
            self.output = self.filtered_dataset['eyey']
            self.z1_num_keys = ["hyhy", "eyhx", "hyhx", "eyhy"]  # Zyx
            self.z2_num_keys = ["hxhx", "eyhy", "hxhy", "eyhx"]  # Zyy
        if self.channel == 'hz':
            self.output = self.filtered_dataset['hzhz']
            self.z1_num_keys = ["hyhy", "hzhx", "hyhx", "hzhy"]  # Zyx
            self.z2_num_keys = ["hxhx", "hzhy", "hxhy", "hzhx"]  # Zyy

    def get_huber_weights(self) -> None:
        """
        Prepares huber weights
        """
        huber_matrix1 = (self.residuals <= self.kmx) * 1
        huber_matrix2 = (self.residuals > self.kmx) * 1
        huber_matrix2 = huber_matrix2 * (self.kmx / self.residuals)
        self.huber_weights = (huber_matrix1 + huber_matrix2).values

    def get_variance(self) -> None:
        """
        Calculates variances.
        Dr. Manoj C. Nair helped with this computation.
        """
        if self.channel == 'ex':
            outout = np.sum(self.filtered_dataset['exex'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhx = np.sum(self.filtered_dataset['exhx'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhy = np.sum(self.filtered_dataset['exhy'].values * self.huber_weights) / np.sum(self.huber_weights)
        if self.channel == 'ey':
            outout = np.sum(self.filtered_dataset['eyey'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhx = np.sum(self.filtered_dataset['eyhx'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhy = np.sum(self.filtered_dataset['eyhy'].values * self.huber_weights) / np.sum(self.huber_weights)
        if self.channel == 'hz':
            outout = np.sum(self.filtered_dataset['hzhz'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhx = np.sum(self.filtered_dataset['hzhx'].values * self.huber_weights) / np.sum(self.huber_weights)
            outhy = np.sum(self.filtered_dataset['hzhy'].values * self.huber_weights) / np.sum(self.huber_weights)

        hxhx = np.sum(self.filtered_dataset['hxhx'].values * self.huber_weights) / np.sum(self.huber_weights)
        hxhy = np.sum(self.filtered_dataset['hxhy'].values * self.huber_weights) / np.sum(self.huber_weights)
        hyhx = np.sum(self.filtered_dataset['hyhx'].values * self.huber_weights) / np.sum(self.huber_weights)
        hyhy = np.sum(self.filtered_dataset['hyhy'].values * self.huber_weights) / np.sum(self.huber_weights)

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

        self.z1_var = np.sqrt(np.real(da))
        self.z2_var = np.sqrt(np.real(db))
        #
        self.predicted_coherency = np.sqrt(coh)
