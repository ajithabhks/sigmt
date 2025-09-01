"""
Utility functions for the metronix
"""

import datetime
import fnmatch
import glob
import os
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd
import xarray as xr
from PyQt5.QtWidgets import QMessageBox

from sigmt.utils.metronix import cal_from_xml


def loadsites(project_dir: str) -> list:
    """
    It lists all folders in a time_series folder in the project, which satisfies two conditions.
    1. Folder should not be empty
    2. Atleast one meas sub folder.
    This is to ensure, it reads proper Metronix sites.

    :param project_dir: Path to the project directory.
    :type project_dir: str
    :return: List of sites
    :rtype: List[str]

    """
    ts_path = os.path.join(project_dir, 'time_series')
    sites = []
    for folder_name in os.listdir(ts_path):
        folder_path = os.path.join(ts_path, folder_name)
        # Check if it is a directory
        if os.path.isdir(folder_path):
            # Check if the directory is not empty
            if os.listdir(folder_path):
                # Check if there is at least one folder starting with 'meas'
                if any(fnmatch.fnmatch(subfolder, 'meas*') for subfolder in os.listdir(folder_path) if
                       os.path.isdir(os.path.join(folder_path, subfolder))):
                    sites.append(folder_name)
    return sites


def list_meas_folders(directory: str) -> list:
    """
    List all folders that start with 'meas' in the given directory.

    This function takes the path to a directory and returns a list of all
    subdirectories that start with 'meas'.

    :param directory: Path to the directory to search.
    :type directory: str
    :return: List of paths to subdirectories starting with 'meas'.
    :rtype: List[str]

    """
    # List to store the values
    meas_folders = []
    # Iterate through the items in the given directory
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        # Check if the item is a directory and starts with "meas"
        if os.path.isdir(item_path) and item.startswith("meas"):
            meas_folders.append(item)
    return meas_folders


def get_sampling_frequency_from_xml(meas_paths: list) -> tuple:
    """
    Get a list of meas paths, loops through the meas folders, reads xml file for sampling frequency
    and chopper status.
    Zips the values, converts to list and returns.

    :param meas_paths: List containing measurement paths.
    :type meas_paths: List[str]
    :return: tuple of lists containing chopper statuses.
    :rtype: Tuple of list[str]

    """
    sampfreq = []
    chopper_value = []
    for meas in meas_paths:
        try:
            allfiles = os.listdir(meas)
            xml_file = next((file for file in allfiles if file.endswith('.xml')), None)
            if xml_file is None:
                raise FileNotFoundError("No XML file in meas folder! Cannot proceed further.")
            xml_file_path = os.path.join(meas, xml_file)
            mytree = ET.parse(xml_file_path)
            root = mytree.getroot()
            sample_freq = root.findall(".//sample_freq")
            sample_freq = float(sample_freq[-1].text)
            sampfreq.append(sample_freq)
            hchopper_value = int(root.find(".//hchopper").text)
            if hchopper_value == 1:
                chopper_status = "chopper_on"
            elif hchopper_value == 0:
                chopper_status = "chopper_off"
            else:
                chopper_status = "unknown"
            chopper_value.append(chopper_status)
        except Exception as e:
            QMessageBox.critical(None, "Error", str(e))
    return sampfreq, chopper_value


def read_ts(meas_path: str, project_setup: dict) -> tuple:
    """
    Reads all time series files (ats) from the measurement path and returns header and time series data as dictionaries.

    :param meas_path: Measurement path
    :type meas_path: str
    :param project_setup: Dictionary containing project information
    :type project_setup: dict
    :return: Header information of time series and Time series data as dictionaries.
    :rtype: Tuple

    """
    header = {}
    ts = {}
    ats_files = glob.glob(os.path.join(meas_path, '*.ats'))
    if not ats_files:
        QMessageBox.information(None, "Information",
                                f"No time series files at {meas_path}. Skipping to next measurement.")
        return {}, {}
    for tspath in ats_files:
        [header_info, ts_data] = read_ats(tspath)
        # Gets start time from header
        start_time = datetime.datetime.utcfromtimestamp(header_info['start'][0])
        # Creates time index
        time_index = pd.date_range(start=start_time, periods=header_info['nsamples'],
                                   freq=str(int(1 / header_info['sfreq'][0] * 1e9)) + 'ns')
        duration = time_index.max() - time_index.min()
        header_info['duration'] = duration  # Saves duration just for information
        header[header_info['ch_type']] = header_info
        # Creates xarray
        if project_setup['processing_mode'] == 'MT + Tipper':
            if header_info['ch_type'] in ['ex', 'ey', 'hx', 'hy', 'hz']:
                ts[header_info['ch_type']] = xr.DataArray(ts_data, coords={'time': time_index}, dims=['time'])
        elif project_setup['processing_mode'] == 'MT Only':
            if header_info['ch_type'] in ['ex', 'ey', 'hx', 'hy']:
                ts[header_info['ch_type']] = xr.DataArray(ts_data, coords={'time': time_index}, dims=['time'])
    return header, ts


def read_ats(filename: str) -> tuple:
    """
    Reads Metronix time series (ATS) data file.
    This function reads both the header information and the time series data from a Metronix ATS file.

    :param filename: List containing measurement paths.
    :type filename: Str
    :return: Header information (Dict) from the ATS file and Time series data (ndarray[Float])
    :rtype: tuple
    """
    with open(filename, 'rb') as f:
        header = {}
        header['length'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()  # Header length
        header['ver'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()  # Header version
        header['nsamples'] = int(np.fromfile(f, dtype=np.int32, count=1)[0])  # Number of samples
        header['sfreq'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()  # sampling frequency, Hz
        # Start time, seconds since 1970
        header['start'] = np.fromfile(f, dtype=np.int32, count=1).tolist()
        header['lsb'] = np.fromfile(
            f, dtype=np.double, count=1).tolist()  # LSB-Value
        header['iGMTOffset'] = np.fromfile(f, dtype=np.int32, count=1).tolist()
        header['rOrigSampleFreq'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['adu_serial_number'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()  # ADU serial number
        header['adu_ADB'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()  # ADU serial number
        header['ch_no'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()  # Channel number
        header['bychopper'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()  # Chopper
        header['ch_type'] = np.fromfile(
            f, dtype=np.int8, count=2).tolist()  # channel type (Hx,Hy,...)
        header['ch_type'] = ''.join([chr(item) for item in header['ch_type']]).lower()
        header['sensor'] = np.fromfile(f, dtype=np.int8, count=6).tolist()
        header['sensor'] = ''.join([chr(item) for item in header['sensor']])
        header['sensor_no'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()  # Sensor serial number
        # x1 coordinate of 1. Dipole (m) - NS
        header['x1'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        # y1 coordinate of 1. Dipole (m) - EW
        header['y1'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        # z1 coordinate of 1. Dipole (m)
        header['z1'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        # x2 coordinate of 1. Dipole (m) - NS
        header['x2'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        # y2 coordinate of 1. Dipole (m) - EW
        header['y2'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        # z2 coordinate of 1. Dipole (m)
        header['z2'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        header['dipole_length'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['dipole_angle'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['rProbeRes'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['rDCOffset'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['rPreGain'] = np.fromfile(f, dtype=np.float32, count=1).tolist()
        header['rPostGain'] = np.fromfile(
            f, dtype=np.float32, count=1).tolist()
        header['lat'] = np.fromfile(f, dtype=np.int32, count=1).tolist()
        header['lon'] = np.fromfile(f, dtype=np.int32, count=1).tolist()
        header['elev'] = np.fromfile(f, dtype=np.int32, count=1).tolist()
        header['byLatLongType'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()
        header['byLatLongType'] = ''.join(
            [chr(item) for item in header['byLatLongType']])
        header['byAddCoordType'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()
        header['byAddCoordType'] = ''.join(
            [chr(item) for item in header['byAddCoordType']])
        header['siRefMedian'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()
        header['dblNorthing'] = np.fromfile(
            f, dtype=np.float64, count=1).tolist()
        header['dblEasting'] = np.fromfile(
            f, dtype=np.float64, count=1).tolist()
        header['byGPSStat'] = np.fromfile(f, dtype=np.int8, count=1).tolist()
        header['byGPSStat'] = ''.join([chr(item)
                                       for item in header['byGPSStat']])
        header['byGPSAccuracy'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()
        header['byGPSAccuracy'] = ''.join(
            [chr(item) for item in header['byGPSAccuracy']])
        header['iUTCOffset'] = np.fromfile(f, dtype=np.int16, count=1).tolist()
        header['achSystemType'] = np.fromfile(
            f, dtype=np.int8, count=12).tolist()
        header['achSystemType'] = ''.join(
            [chr(item) for item in header['achSystemType']])
        header['achSurveyHeaderName'] = np.fromfile(
            f, dtype=np.int8, count=12).tolist()
        header['achSurveyHeaderName'] = ''.join(
            [chr(item) for item in header['achSurveyHeaderName']])
        header['achMeasType'] = np.fromfile(f, dtype=np.int8, count=4).tolist()
        header['achMeasType'] = ''.join(
            [chr(item) for item in header['achMeasType']])
        header['DCOffsetCorrValue'] = np.fromfile(
            f, dtype=np.float64, count=1).tolist()
        header['DCOffsetCorrOn'] = np.fromfile(
            f, dtype=np.int8, count=1).tolist()
        header['InputDivOn'] = np.fromfile(f, dtype=np.int8, count=1).tolist()
        header['bit_indicator'] = np.fromfile(
            f, dtype=np.int16, count=1).tolist()
        header['achSelfTestResult'] = np.fromfile(
            f, dtype=np.int8, count=2).tolist()
        header['achSelfTestResult'] = ''.join(
            [chr(item) for item in header['achSelfTestResult']])
        header['numslice'] = np.fromfile(f, dtype=np.uint16, count=1).tolist()
        header['siCalFreqs'] = np.fromfile(
            f, dtype=np.uint16, count=1).tolist()
        f.seek(header['length'][0])
        ts = np.fromfile(f, dtype=np.int32, count=header['nsamples'])
        ts = ts * header['lsb'][0]
        return header, ts


def read_calibration_from_xml(meas_path: str) -> dict:
    """
    Reads calibration data from xml files.

    :param meas_path: Path to the measurement file, for eg: 'D:/project\\time_series\\KL33A\\meas_2019-08-29_14-30-00'
    :type meas_path: str
    :return: The Dictionary contains calibration data, ['coil_number'], within that 'ChoppOn' and 'ChoppOff'
    :rtype: dict

    """
    calibration_data = {}
    xml_files = []
    for file in os.listdir(meas_path):
        if file.endswith(".xml"):
            xml_files.append(file)
    if len(xml_files) == 0:
        print("No XML found.")
        calibration_data = None
    else:
        print("XML found. Attempting to read calibration data.")
        xml_file = os.path.join(meas_path, xml_files[0])
        channels_info = cal_from_xml.extract_channel_info(xml_file)
        # Filter out channels with types 'Ex' and 'Ey'
        filtered_channels = cal_from_xml.filter_channels(channels_info, ['Ex', 'Ey'])
        filtered_channel_ids = [channel[0] for channel in filtered_channels]
        for filtered_channel_id, val in zip(filtered_channel_ids, filtered_channels):
            calibration_data[val[2]] = cal_from_xml.read_cal_data(xml_file, filtered_channel_id)
        print("Calibration data reading from XML done!")

    return calibration_data


def prepare_ts_from_h5(h5file_path: str, ts_key: str) -> dict:
    """
    Read time series from the h5 file.

    :param h5file_path: Path to the database file
    :type h5file_path: str
    :param ts_key: Key value for to access from the h5 group. Eg: 'ts_0'
    :type ts_key: str
    :return: Dictionary containing time series for 'Ex', 'Ey', ...
    :rtype: dict

    """
    ts = {}
    with h5py.File(h5file_path, 'r') as f:
        channel_keys = f[ts_key].keys()
        for channel in channel_keys:
            ts[channel] = f[ts_key][channel][:]

    return ts
