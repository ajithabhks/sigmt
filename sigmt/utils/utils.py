"""
Place to store utility functions.
"""
from typing import Dict, Any

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
    with open(file_path, 'r', encoding='utf-8') as yaml_file:
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

    for _, sub_dict in header.items():
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
    term = nofsamples / (20 * i)
    ffts.append(cfft)
    while cfft <= term:
        i = i + 1
        cfft = 256 * (2 ** i)
        term = nofsamples / (20 * i)
        if cfft <= term:
            ffts.append(cfft)

    window_length = min(ffts[-1], 65536)
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


def targetfreq(fs, cr, fftlength, periods_per_decade):
    """
    It returns target frequencies corresponding to sampling frequency.

    :param fs: Sampling frequency of measurement
    :type fs: float
    :param cr: Parzen window radius
    :type cr: float
    :param fftlength: FFT length
    :type fftlength: int
    :param freq_per_decade: Frequencies per decade required
    :type freq_per_decade: int

    :returns: A numpy array (1D) of float (shape: n,) which is a list of target frequencies.
    :rtype: np.ndarray

    """
    start_period = -5
    stop_period = 5
    periods_per_decade = 12
    ftable = np.logspace(start_period, stop_period, int(
        (stop_period - start_period) * periods_per_decade + 1))
    ftable = 1 / ftable

    fr = cr * ftable  # bandwidth of parzen window - oneside from ft
    totalbandwidth = fr * 2  # bandwidth of parzen window - two side from ft
    # nof spectra in parzen window for each ft
    dof = totalbandwidth / (fs / fftlength)

    maximum = fs / 2
    maximum = min(maximum, 15000)
    fmax = max(ftable[ftable < maximum])  # ft_max, following nyquist

    fmaxindex = np.where(ftable == fmax)[0][0]  # index in ftable
    min_required = dof[fmaxindex] * .01  # min dof required for ft_min

    dofmin = min(dof[dof > min_required])  # finding in dof
    dofminindex = np.where(dof == dofmin)[0][0]  # finding index

    # fmin = ftable[dofminindex] #ft_min

    while dof[dofminindex] < 10:  # Making sure, atleast 10 spectral lines are averaged
        dofminindex = dofminindex - 1

    # fmin = ftable[dofminindex] #Fixing ft_min

    ftlist = ftable[fmaxindex:dofminindex + 1]  # Making ftlist

    return ftlist
