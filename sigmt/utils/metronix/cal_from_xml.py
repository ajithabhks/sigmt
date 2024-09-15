"""
Functions to read calibration data from the xml files.
"""

import xml.etree.ElementTree as ET

import numpy as np


def extract_channel_info(xmlfile: str) -> list:
    """
    Extracts number of channels in the xml file. It gives the channel IDs. It is required to identify
    coil if some calibration sections are missing from xml.

    :param xmlfile: Path to the xml file
    :type xmlfile: str

    :return: Contains a list of tuples. Eg tuple: ('7', 'Hx', '300'). 7: channel id, 300: sensor number.
    :rtype: list
    """
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    channels_info = []

    for configuration in root.iter('configuration'):
        # Iterate over each <channel> within the <configuration>
        for channel in configuration.findall('channel'):
            channel_id = channel.get('id')
            channel_type = channel.find('channel_type').text if channel.find('channel_type') is not None else None
            sensor_sernum = channel.find('sensor_sernum').text if channel.find('sensor_sernum') is not None else None
            channels_info.append((channel_id, channel_type, sensor_sernum))

    return channels_info


def filter_channels(data: list, excluded_types: list) -> list:
    """
    It is to exclude some channels which are not required. For eg: 'Ex', 'Ey'

    :param data: It is the list of the tuples. Eg: ('5', 'Ex', '0')
    :type data: list
    :param excluded_types: List of strings, to be excluded. Eg: 'Ex', 'Ey'
    :type excluded_types: list

    :return: List of channels (data above) after exclusion.
    :rtype: list
    """
    return [channel for channel in data if channel[1] not in excluded_types]


def read_cal_data(xmlfile: str, channel_id: str) -> dict:
    """
    Reads freq, amplitude, phase for chopper on and off, for a coil in xml file.

    :param xmlfile: Path to the xml file
    :type xmlfile: str
    :param channel_id: Channel id
    :type channel_id: str

    :return: Dictionary containing 'ChoppOn' and 'ChoppOff' data.
    :rtype: dict
    """
    mytree = ET.parse(xmlfile)
    root = mytree.getroot()
    #

    calsensors = root[2]

    calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='on']/c1")

    caldata = {}
    freq = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='on']/c1"):
        freq.append(float(x.text))

    mag = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='on']/c2"):
        mag.append(float(x.text))

    phase = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='on']/c3"):
        phase.append(float(x.text))

    freq = np.asarray(freq).reshape(-1, 1)
    mag = np.asarray(mag).reshape(-1, 1)
    phase = np.asarray(phase).reshape(-1, 1)
    caldata['ChoppOn'] = np.concatenate((freq, mag, phase), axis=1)

    freq = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='off']/c1"):
        freq.append(float(x.text))

    mag = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='off']/c2"):
        mag.append(float(x.text))

    phase = []
    for x in calsensors.findall(f".//channel[@id='{channel_id}']//caldata[@chopper='off']/c3"):
        phase.append(float(x.text))

    freq = np.asarray(freq).reshape(-1, 1)
    mag = np.asarray(mag).reshape(-1, 1)
    phase = np.asarray(phase).reshape(-1, 1)
    caldata['ChoppOff'] = np.concatenate((freq, mag, phase), axis=1)

    return caldata
