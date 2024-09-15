"""
Functions to plot output curves
"""

import matplotlib
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

plt.rcParams['figure.max_open_warning'] = 50


def plot_mt_app_res(dataset: xr.Dataset) -> None:
    """
    Plots apparent resistivity and phase

    :param dataset: Dataset containing estimates
    :type dataset: xr.Dataset
    :return: None
    :rtype: None

    """
    if dataset:
        # Apparant resistivities and phase
        ftlist = dataset['frequency'].values
        zxy = dataset['zxy'].values
        zyx = dataset['zyx'].values
        zxy_var = dataset['zxy_var'].values
        zyx_var = dataset['zyx_var'].values

        rho_xy = (0.2 / ftlist) * (abs(zxy) ** 2)
        rho_yx = (0.2 / ftlist) * (abs(zyx) ** 2)
        phase_xy = np.degrees(np.arctan(zxy.imag / zxy.real))
        phase_yx = np.degrees(np.arctan(zyx.imag / zyx.real))
        #
        #
        # Errors for app. resistivity and phase
        err_rxy = (0.4 / ftlist) * abs(zxy) * zxy_var
        err_ryx = (0.4 / ftlist) * abs(zyx) * zyx_var
        err_pxy = np.degrees(zxy_var / abs(zxy))
        err_pyx = np.degrees(zyx_var / abs(zyx))
        #
        # Plotting section
        # Plot App. res & Phase
        plt.figure()
        plt.subplot(211)
        plt.scatter(ftlist, rho_xy, c='r', s=10, label='XY')
        plt.scatter(ftlist, rho_yx, c='b', s=10, label='YX')
        plt.errorbar(ftlist, rho_xy, yerr=err_rxy, ecolor='r', fmt="none")
        plt.errorbar(ftlist, rho_yx, yerr=err_ryx, ecolor='b', fmt="none")
        plt.xscale('log')
        plt.yscale('log')
        if max(ftlist) < 15000 and min(ftlist) > 0.001:
            plt.xlim((15000, 0.001))
        else:
            plt.xlim((max(ftlist) + 10, min(ftlist) - 10))
            if min(ftlist) > 0.001:
                plt.xlim((max(ftlist) + 10, 0.001))
        ymax_magnitude = int(np.ceil(np.log10(max(max(rho_xy), max(rho_yx)))))
        ymin_magnitude = int(np.floor(np.log10(min(min(rho_xy), min(rho_yx)))))
        plt.ylim(10 ** (ymin_magnitude - 2), 10 ** (ymax_magnitude + 2))
        plt.yticks([10 ** i for i in range(ymin_magnitude - 2, ymax_magnitude + 3)])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('App. Res. (Ohm.m.)')
        plt.legend()
        plt.grid(which='both', linestyle='-.', linewidth=0.4)
        plt.subplot(212)
        plt.scatter(ftlist, phase_xy, c='r', s=10)
        plt.scatter(ftlist, phase_yx, c='b', s=10)
        plt.errorbar(ftlist, phase_xy, yerr=err_pxy, ecolor='r', fmt="none")
        plt.errorbar(ftlist, phase_yx, yerr=err_pyx, ecolor='b', fmt="none")
        plt.xscale('log')
        if max(ftlist) < 15000 and min(ftlist) > 0.001:
            plt.xlim((15000, 0.001))
        else:
            plt.xlim((max(ftlist) + 10, min(ftlist) - 10))
            if min(ftlist) > 0.001:
                plt.xlim((max(ftlist) + 10, 0.001))
        if (min(phase_xy) > 0 and min(phase_yx) > 0) and (max(phase_xy) < 90 and max(phase_yx) < 90):
            plt.ylim((0, 90))
            plt.yticks([0, 15, 30, 45, 60, 75, 90])
        else:
            ylim1 = min(min(phase_xy), min(phase_yx))
            ylim2 = max(max(phase_xy), max(phase_yx))
            plt.ylim(ylim1 - 10, ylim2 + 10)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Phase (Deg.)')
        plt.grid(which='both', linestyle='-.', linewidth=0.4)
        plt.show()
    #


def plot_tipper(dataset: xr.Dataset) -> None:
    """
    Plots tipper

    :param dataset: Dataset containing estimates
    :type dataset: xr.Dataset
    :return: None
    :rtype: None

    """
    if dataset:
        ftlist = dataset['frequency'].values
        tzx = dataset['tzx'].values
        txy = dataset['tzy'].values
        #
        tx_a = np.sqrt((np.real(tzx) ** 2) + (np.imag(tzx) ** 2))
        ty_a = np.sqrt((np.real(txy) ** 2) + (np.imag(txy) ** 2))
        tx_p = np.degrees(np.arctan2(tzx.imag, tzx.real))
        ty_p = np.degrees(np.arctan2(txy.imag, txy.real))
        # Plotting section
        plt.figure()
        plt.subplot(211)
        plt.scatter(ftlist, tx_a, c='r', s=10, label='Tx')
        plt.scatter(ftlist, ty_a, c='b', s=10, label='Ty')
        plt.ylim(0, 1)
        plt.xscale('log')
        if max(ftlist) < 15000 and min(ftlist) > 0.001:
            plt.xlim((15000, 0.001))
        else:
            plt.xlim((max(ftlist) + 10, min(ftlist) - 10))
            if min(ftlist) > 0.001:
                plt.xlim((max(ftlist) + 10, 0.001))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Tipper Amplitude')
        plt.grid(which='both', linestyle='-.', linewidth=0.4)
        plt.legend()
        plt.subplot(212)
        plt.scatter(ftlist, tx_p, c='r', s=10)
        plt.scatter(ftlist, ty_p, c='b', s=10)
        plt.xscale('log')
        if max(ftlist) < 15000 and min(ftlist) > 0.001:
            plt.xlim((15000, 0.001))
        else:
            plt.xlim((max(ftlist) + 10, min(ftlist) - 10))
            if min(ftlist) > 0.001:
                plt.xlim((max(ftlist) + 10, 0.001))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Tipper Phase')
        plt.grid(which='both', linestyle='-.', linewidth=0.4)
        plt.show()


def plot_coherency(dataset: xr.Dataset) -> None:
    """
    Plots coherencies

    :param dataset: Dataset containing estimates
    :type dataset: xr.Dataset
    :return: None
    :rtype: None

    """
    if dataset:
        plt.figure()
        if 'coh_ex' in dataset:
            plt.scatter(dataset['frequency'].values, dataset['coh_ex'].values, c='r', label='Ex')
        if 'coh_ey' in dataset:
            plt.scatter(dataset['frequency'].values, dataset['coh_ey'].values, c='b', label='Ey')
        if 'coh_hz' in dataset:
            plt.scatter(dataset['frequency'].values, dataset['coh_hz'].values, c='g', label='Hz')
        plt.xscale('log')
        plt.xlim((10000, 0.001))
        plt.ylim(0, 1)
        plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Predicted Coherency')
        plt.grid(which='both', linestyle='-.', linewidth=0.4)
        plt.legend()
        plt.show()


def plot_coherencies_all(dataset: xr.Dataset) -> None:
    """
    Plots coherency values for all frequencies

    :param dataset: Dataset containing band averaged values
    :type dataset: xr.Dataset
    :return: None
    :rtype: None

    """
    if dataset:
        cdict = {'red': ((0.0, 0.0, 0.0),
                         (0.1, 0.5, 0.5),
                         (0.2, 0.0, 0.0),
                         (0.4, 0.2, 0.2),
                         (0.6, 0.0, 0.0),
                         (0.8, 1.0, 1.0),
                         (1.0, 1.0, 1.0)),
                 'green': ((0.0, 0.0, 0.0),
                           (0.1, 0.0, 0.0),
                           (0.2, 0.0, 0.0),
                           (0.4, 1.0, 1.0),
                           (0.6, 1.0, 1.0),
                           (0.8, 1.0, 1.0),
                           (1.0, 0.0, 0.0)),
                 'blue': ((0.0, 0.0, 0.0),
                          (0.1, 0.5, 0.5),
                          (0.2, 1.0, 1.0),
                          (0.4, 1.0, 1.0),
                          (0.6, 0.0, 0.0),
                          (0.8, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}
        my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

        num_coh = 0
        coh_keys = []
        if 'coh_ex' in dataset:
            num_coh = num_coh + 1
            coh_keys.append('coh_ex')
        if 'coh_ey' in dataset:
            num_coh = num_coh + 1
            coh_keys.append('coh_ey')
        if 'coh_hz' in dataset:
            num_coh = num_coh + 1
            coh_keys.append('coh_hz')
        #
        for freq in dataset['frequency']:
            # Select data for the current frequency
            data_for_freq = dataset.sel(frequency=freq)
            if num_coh == 1:
                if coh_keys[0] == 'coh_ex':
                    x = 'zxy_single'
                    y = 'zxy_single'
                if coh_keys[0] == 'coh_ey':
                    x = 'zyx_single'
                    y = 'zyx_single'
                if coh_keys[0] == 'coh_hz':
                    x = 'tzx_single'
                    y = 'tzx_single'
                plt.figure()
                plt.scatter(data_for_freq[x].real, data_for_freq[y].imag, c=data_for_freq[coh_keys[0]], cmap=my_cmap,
                            vmin=0, vmax=1)
                plt.colorbar(label=coh_keys[0])  # Add color bar for reference
                plt.title(f"Coherency ({coh_keys[0]}) plot for {freq.values} Hz", fontsize=12)
                plt.xlabel(x)
                plt.ylabel(y)
            elif num_coh == 2:
                fig, axes = plt.subplots(1, 2, figsize=(10, 5))
                for i in range(2):
                    if coh_keys[i] == 'coh_ex':
                        x = 'zxy_single'
                        y = 'zxy_single'
                    if coh_keys[i] == 'coh_ey':
                        x = 'zyx_single'
                        y = 'zyx_single'
                    if coh_keys[i] == 'coh_hz':
                        x = 'tzx_single'
                        y = 'tzx_single'
                    sc = axes[i].scatter(data_for_freq[x].real, data_for_freq[y].imag, c=data_for_freq[coh_keys[i]],
                                         cmap=my_cmap, vmin=0, vmax=1)
                    fig.colorbar(sc, ax=axes[i], label=coh_keys[i])
                    axes[i].set_title(coh_keys[i])
                    axes[i].set_xlabel(x)
                    axes[i].set_ylabel(y)
                fig.suptitle(f"Coherency plots for {freq.values} Hz", fontsize=12)
            elif num_coh == 3:
                fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 3 horizontal subplots
                for i in range(3):
                    if coh_keys[i] == 'coh_ex':
                        x = 'zxy_single'
                        y = 'zxy_single'
                    if coh_keys[i] == 'coh_ey':
                        x = 'zyx_single'
                        y = 'zyx_single'
                    if coh_keys[i] == 'coh_hz':
                        x = 'tzx_single'
                        y = 'tzx_single'
                    sc = axes[i].scatter(data_for_freq[x].real, data_for_freq[y].imag, c=data_for_freq[coh_keys[i]],
                                         cmap=my_cmap, vmin=0, vmax=1)
                    fig.colorbar(sc, ax=axes[i], label=coh_keys[i])
                    axes[i].set_title(coh_keys[i])
                    axes[i].set_xlabel(x)
                    axes[i].set_ylabel(y)
                fig.suptitle(f"Coherency plots for {freq.values} Hz", fontsize=12)

            plt.tight_layout()
            plt.show()


def plot_pd_all(dataset: xr.Dataset) -> None:
    """
    Plots polarization directions values for all frequencies

    :param dataset: Dataset containing band averaged values
    :type dataset: xr.Dataset
    :return: None
    :rtype: None

    """
    if dataset:
        cdict = {'red': ((0.0, 0.0, 0.0),
                         (0.1, 0.5, 0.5),
                         (0.2, 0.0, 0.0),
                         (0.4, 0.2, 0.2),
                         (0.6, 0.0, 0.0),
                         (0.8, 1.0, 1.0),
                         (1.0, 1.0, 1.0)),
                 'green': ((0.0, 0.0, 0.0),
                           (0.1, 0.0, 0.0),
                           (0.2, 0.0, 0.0),
                           (0.4, 1.0, 1.0),
                           (0.6, 1.0, 1.0),
                           (0.8, 1.0, 1.0),
                           (1.0, 0.0, 0.0)),
                 'blue': ((0.0, 0.0, 0.0),
                          (0.1, 0.5, 0.5),
                          (0.2, 1.0, 1.0),
                          (0.4, 1.0, 1.0),
                          (0.6, 0.0, 0.0),
                          (0.8, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}
        my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)
        afont = {'fontname': 'Arial'}
        yticks = np.arange(-90, 100, 20)

        for freq in dataset['frequency']:
            data_for_freq = dataset.sel(frequency=freq)
            time_array = np.arange(data_for_freq.coords['time_window'].size)
            plt.figure()
            plt.subplot(211)
            plt.suptitle(f"Polarization directions for {freq.values} Hz", fontsize=12)
            plt.scatter(time_array, data_for_freq['alpha_h'].values)
            plt.ylim(-90, 90)
            plt.ylabel(r'$\alpha_H$ (deg.)')
            plt.xticks(**afont, fontsize=12)
            plt.yticks(yticks, **afont, fontsize=12)
            plt.subplot(212)
            plt.scatter(time_array, data_for_freq['alpha_e'].values)
            plt.ylim(-90, 90)
            plt.ylabel(r'$\alpha_E$ (deg.)')
            plt.xlabel('Time window')
            plt.xticks(**afont, fontsize=12)
            plt.yticks(yticks, **afont, fontsize=12)
            plt.show(block=False)
