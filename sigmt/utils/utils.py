"""
Place to store utility functions.
"""
from typing import Dict, Any, Optional

import numpy as np
import yaml


def read_yaml_file(file_path: str) -> Dict[str, Any]:
    """
    Reads the yaml file and returns as a dictionary.

    :param file_path: Path to yaml file
    :type file_path: str
    :return: Content in the yaml file
    :rtype: dict
    """
    with open(file_path, 'r') as yaml_file:
        data = yaml.load(yaml_file, Loader=yaml.FullLoader)
    return data


def get_nsamples(header: dict) -> int:
    """
    Takes header dictionary, gets all nsamples.
    Finds the mostly occurring number and returns.
    
    :param header: Header information of all loaded time series.
    :type header: dict
    :return: Most common number in nsamples
    """
    nsamples_list = []

    for key, sub_dict in header.items():
        nsamples = sub_dict.get('nsamples')
        if nsamples is not None:
            nsamples_list.append(nsamples)
        else:
            # If 'nsamples' key is not found in the sub-dictionary, check its values
            for sub_key, value in sub_dict.items():
                if sub_key not in ['Rx', 'Ry'] and isinstance(value, dict):
                    nested_nsamples = value.get('nsamples')
                    if nested_nsamples is not None:
                        nsamples_list.append(nested_nsamples)

    most_common_number = max(set(nsamples_list), key=nsamples_list.count)
    return most_common_number


def get_fftlength(nofsamples: int) -> int:
    """
    Returns a suggested FFT length based on the no. of samples.
    :param nofsamples: No. of samples
    :type nofsamples: int
    :return: Suggested FFT length
    :rtype: int
    """
    # Based on Borah & Patro, 2015
    i = 1
    ffts = [256]
    cfft = 256 * (2 ** i)
    term = nofsamples / (20 * i);
    ffts.append(cfft)
    while cfft <= term:
        i = i + 1
        cfft = 256 * (2 ** i)
        term = nofsamples / (20 * i)
        if cfft <= term:
            ffts.append(cfft)
    window_length = ffts[-1]
    if window_length > 65536:
        window_length = 65536
    return window_length


def get_parzen(fs: float) -> float:
    """
    Returns a suggested parzen window radius based on the sampling frequency.

    :param fs: Sampling frequency
    :type fs: float
    :return: Parzen window radius
    :rtype: float
    """
    if fs >= 512:
        parzen_radius = 0.25
    elif fs > 32:
        parzen_radius = 0.4
    else:  # fs <= 32
        parzen_radius = 0.9
    return parzen_radius


def reshape_array_with_overlap(window_length: int, overlap: int, data: np.ndarray) -> np.ndarray:
    """
    To reshape a numpy array into different windows with 50% overlap.

    :param window_length: Window length that arrays to be divided.
    :type window_length: int
    :param overlap: Overlap required for the reshaping.
    :type overlap: int
    :param data: The data array (1D)
    :type data: np.ndarray

    :return: 2D numpy array, rows= no. of windows, columns = window length
    :rtype: np.ndarray
    """
    # Parameters
    overlap = int(window_length * (overlap / 100))

    # Calculate the number of windows
    num_windows = (len(data) - overlap) // (window_length - overlap)

    # Initialize the 2D array
    result_array = np.empty((num_windows, window_length))

    # Extract windows
    for i in range(num_windows):
        start_idx = i * (window_length - overlap)
        end_idx = start_idx + window_length
        result_array[i] = data[start_idx:end_idx]

    result_array = result_array.T

    return result_array


def get_target_frequency_list(sampling_frequency: float,
                              parzen_window_radius: float,
                              fft_length: int,
                              periods_per_decade: int,
                              start_period: int = -5,
                              stop_period: int = 5
                              ) -> np.ndarray:
    """
    Computes a list of target frequencies based on the given processing parameters.

    :param sampling_frequency: Sampling frequency of the measurement (Hz).
    :type sampling_frequency: float
    :param parzen_window_radius: Radius of the Parzen window.
    :type parzen_window_radius: float
    :param fft_length: Length of the FFT (Fast Fourier Transform).
    :type fft_length: int
    :param periods_per_decade: Number of period/frequency points per decade.
    :type periods_per_decade: int
    :param start_period: Exponent for the start of the period range (10^start_period).
                         Default is -5, corresponding to 10⁻⁵ Hz.
    :type start_period: int
    :param stop_period: Exponent for the end of the period range (10^stop_period).
                        Default is 5, corresponding to 10⁵ Hz.
    :type stop_period: int

    :returns: A numpy array (1D) of float (shape: n, ) which is a list of target frequencies.
    :rtype: np.ndarray

    """

    ftable = np.logspace(start_period, stop_period, int(
        (stop_period - start_period) * periods_per_decade + 1))
    ftable = 1 / ftable

    fr = parzen_window_radius * ftable  # bandwidth of parzen window - oneside from ft
    total_bandwidth = fr * 2  # bandwidth of parzen window - two side from ft
    # nof spectra in parzen window for each ft
    dof = total_bandwidth / (sampling_frequency / fft_length)

    maximum = sampling_frequency / 2
    if maximum > 15000:
        maximum = 15000
    f_max = max(ftable[ftable < maximum])  # ft_max, following nyquist

    f_max_index = np.where(ftable == f_max)[0][0]  # index in ftable
    min_required = dof[f_max_index] * .01  # min dof required for ft_min

    dof_min = min(dof[dof > min_required])  # finding in dof
    dof_min_index = np.where(dof == dof_min)[0][0]  # finding index

    while dof[dof_min_index] < 10:  # Making sure, atleast 10 spectral lines are averaged
        dof_min_index = dof_min_index - 1

    ft_list = ftable[f_max_index:dof_min_index + 1]  # Making ft_list

    return ft_list
