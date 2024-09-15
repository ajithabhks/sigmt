import numpy as np
import xarray as xr


def cohex(bandavg: xr.Dataset) -> np.ndarray:
    """
    Computes coherency values for ex output channel.

    :param bandavg: Dataset containing the auto- and cross- spectra values and
    impedance values for all time windows at all target frequencies.
    :type bandavg: xr.Dataset
    :return: Returns coherency values for all events.
    :rtype: np.ndarray

    """
    # Ex predicted
    zxx = bandavg['zxx_single']
    zxy = bandavg['zxy_single']
    zpz = zxx * np.conj(bandavg['exhx']) + zxy * np.conj(bandavg['exhy'])
    zpx = zxx * bandavg['hxhx'] + zxy * bandavg['hyhx']
    zpy = zxx * bandavg['hxhy'] + zxy * bandavg['hyhy']
    zpzp = zxx * np.conj(zpx) + zxy * np.conj(zpy)
    zpzp = zpzp * bandavg['exex']

    coh = np.where(np.abs(zpzp) > 0, zpz / np.sqrt(zpzp), 1 + 1j)
    coh = np.abs(coh)

    coh = np.where(coh > 1.0, 1 / coh, coh)
    return coh


def cohey(bandavg: xr.Dataset) -> np.ndarray:
    """
    Computes coherency values for ey output channel.

    :param bandavg: Dataset containing the auto- and cross- spectra values and
    impedance values for all time windows at all target frequencies.
    :type bandavg: xr.Dataset
    :return: Returns coherency values for all events.
    :rtype: np.ndarray

    """
    # Ey predicted
    zyx = bandavg['zyx_single']
    zyy = bandavg['zyy_single']
    zpz = zyx * np.conj(bandavg['eyhx']) + zyy * np.conj(bandavg['eyhy'])
    zpx = zyx * bandavg['hxhx'] + zyy * bandavg['hyhx']
    zpy = zyx * bandavg['hxhy'] + zyy * bandavg['hyhy']
    zpzp = zyx * np.conj(zpx) + zyy * np.conj(zpy)
    zpzp = zpzp * bandavg['eyey']

    coh = np.where(np.abs(zpzp) > 0, zpz / np.sqrt(zpzp), 1 + 1j)
    coh = np.abs(coh)

    coh = np.where(coh > 1.0, 1 / coh, coh)
    return coh


def cohhz(bandavg: xr.Dataset) -> np.ndarray:
    """
    Computes coherency values for hz output channel.

    :param bandavg: Dataset containing the auto- and cross- spectra values and
    impedance values for all time windows at all target frequencies.
    :type bandavg: xr.Dataset
    :return: Returns coherency values for all events.
    :rtype: np.ndarray

    """
    # Ex predicted
    tzx = bandavg['tzx_single']
    tzy = bandavg['tzy_single']
    zpz = tzx * np.conj(bandavg['hzhx']) + tzy * np.conj(bandavg['hzhy'])
    zpx = tzx * bandavg['hxhx'] + tzy * bandavg['hyhx']
    zpy = tzx * bandavg['hxhy'] + tzy * bandavg['hyhy']
    zpzp = tzx * np.conj(zpx) + tzy * np.conj(zpy)
    zpzp = zpzp * bandavg['hzhz']

    coh = np.where(np.abs(zpzp) > 0, zpz / np.sqrt(zpzp), 1 + 1j)
    coh = np.abs(coh)

    coh = np.where(coh > 1.0, 1 / coh, coh)
    return coh


def pdvalues(bandavg: xr.Dataset) -> tuple:
    """
    Computes polarization directions for electric and magnetic field.

    :param bandavg: Dataset containing the auto- and cross- spectra values and
    impedance values for all time windows at all target frequencies.
    :type bandavg: xr.Dataset
    :return: Tuple containing polarization directions for electric field and magnetic field.
    :rtype: tuple

    """
    hxhx = bandavg['hxhx']
    hyhy = bandavg['hyhy']
    hxhy = bandavg['hxhy']
    alpha_h = np.arctan(2 * np.real(hxhy) / (hxhx - hyhy))
    alpha_h_deg = np.degrees(np.real(alpha_h))

    exex = bandavg['exex']
    eyey = bandavg['eyey']
    exey = bandavg['exey']
    alpha_e = np.arctan(2 * np.real(exey) / (exex - eyey))
    alpha_e_deg = np.degrees(np.real(alpha_e))
    return alpha_h_deg, alpha_e_deg
