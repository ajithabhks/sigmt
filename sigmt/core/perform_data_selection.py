import numpy as np
import xarray as xr


def perform_coh_thresh(bandavg_dataset, coh_thresh, min_percent, channels):
    """
    Doc
    """
    for channel in channels:
        min_windows = int(
            np.ceil(len(bandavg_dataset[f'coh_{channel}'].coords['time_window']) * (min_percent / 100)))
        for ft in bandavg_dataset[f'coh_{channel}'].coords['frequency']:
            ct = coh_thresh
            coh = bandavg_dataset[f'coh_{channel}'].sel(frequency=ft)
            while np.sum(coh >= ct) < min_windows:
                ct = ct - 0.01
            selection_array = coh > ct
            bandavg_dataset[f'{channel}_selection_coh'].loc[dict(frequency=ft)] = selection_array

    return bandavg_dataset


def clear_coh_thresh(bandavg_dataset):
    """
    Doc
    """
    if 'ex_selection_coh' in bandavg_dataset:
        bandavg_dataset['ex_selection_coh'] = xr.DataArray(
            np.full(bandavg_dataset['ex_selection_coh'].shape, True),
            coords=bandavg_dataset.coords,
            dims=bandavg_dataset.dims
        )
    if 'ey_selection_coh' in bandavg_dataset:
        bandavg_dataset['ey_selection_coh'] = xr.DataArray(
            np.full(bandavg_dataset['ey_selection_coh'].shape, True),
            coords=bandavg_dataset.coords,
            dims=bandavg_dataset.dims
        )
    if 'hz_selection_coh' in bandavg_dataset:
        bandavg_dataset['hz_selection_coh'] = xr.DataArray(
            np.full(bandavg_dataset['hz_selection_coh'].shape, True),
            coords=bandavg_dataset.coords,
            dims=bandavg_dataset.dims
        )
    return bandavg_dataset


def perform_pd_selection(bandavg_dataset, component, pd_min, pd_max):
    """
    Doc
    """
    if component == 'Electric':
        for ft in bandavg_dataset['alpha_e'].coords['frequency']:
            alpha_e = bandavg_dataset['alpha_e'].sel(frequency=ft)
            selection_array = (alpha_e < pd_min) & (alpha_e > pd_max)
            bandavg_dataset['alpha_e_selection'].loc[dict(frequency=ft)] = selection_array
    elif component == 'Magnetic':
        for ft in bandavg_dataset['alpha_h'].coords['frequency']:
            alpha_h = bandavg_dataset['alpha_h'].sel(frequency=ft)
            selection_array = (alpha_h < pd_min) & (alpha_h > pd_max)
            bandavg_dataset['alpha_h_selection'].loc[dict(frequency=ft)] = selection_array
    return bandavg_dataset


def clear_pd_selection(bandavg_dataset):
    """
    Doc
    """
    if 'alpha_e_selection' in bandavg_dataset:
        bandavg_dataset['alpha_e_selection'] = xr.DataArray(
            np.full(bandavg_dataset['ex'].shape, True),
            coords=bandavg_dataset.coords,
            dims=bandavg_dataset.dims
        )

    if 'alpha_h_selection' in bandavg_dataset:
        bandavg_dataset['alpha_h_selection'] = xr.DataArray(
            np.full(bandavg_dataset['hx'].shape, True),
            coords=bandavg_dataset.coords,
            dims=bandavg_dataset.dims
        )
    return bandavg_dataset
