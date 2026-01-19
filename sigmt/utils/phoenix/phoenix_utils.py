"""
Utility functions for the phoenix
"""

import os
import pathlib
import re
from math import ceil
from typing import List, Set, Tuple, Dict, Optional

import numpy as np
from PhoenixGeoPy.Reader import TimeSeries as PhoenixReader

from sigmt.utils import utils
from sigmt.utils.phoenix import phoenix_readers


def load_sites(project_dir: str) -> List[str]:
    """
    List valid Phoenix sites found under the `time_series` directory.

    A site is considered valid if:
    1. It is a directory
    2. It is not empty
    3. It contains at least one recording folder matching the pattern:
       <digits>_YYYY-MM-DD-HHMMSS

    :param project_dir: Path to the project directory
    :return: List of valid site names
    """

    recording_pattern = re.compile(
        r"^\d+_\d{4}-\d{2}-\d{2}-\d{6}$"
    )

    ts_path = os.path.join(project_dir, "time_series")

    if not os.path.isdir(ts_path):
        raise FileNotFoundError(f"'time_series' directory not found: {ts_path}")

    sites: List[str] = []

    for site_name in os.listdir(ts_path):
        site_path = os.path.join(ts_path, site_name)

        if not os.path.isdir(site_path):
            continue

        try:
            entries = os.listdir(site_path)
        except PermissionError:
            continue

        if not entries:
            continue

        # Check for at least one valid recording folder
        has_recording = any(
            os.path.isdir(os.path.join(site_path, entry))
            and recording_pattern.match(entry)
            for entry in entries
        )

        if has_recording:
            sites.append(site_name)

    return sorted(sites)


def get_sampling_rate_list(recording_path):
    """
    Return a sorted list of unique sampling rates found in a recording folder.

    Mapping:
        td_24k  -> 24000
        td_2400 -> 2400
        td_150  -> 150
        td_30   -> 30
    """
    extension_to_sampling_rate = {
        "td_24k": 24000,
        "td_2400": 2400,
        "td_150": 150,
        "td_30": 30,
    }

    file_extensions = list_unique_td_extensions(recording_path=recording_path)

    sampling_rates = {
        str(extension_to_sampling_rate[ext])
        for ext in file_extensions
        if ext in extension_to_sampling_rate
    }

    return sorted(sampling_rates)


def sampling_rate_to_extension(sampling_rate):
    """
    Convert a sampling rate to its corresponding td_* extension.

    Examples:
        24000 -> "td_24k"
        2400  -> "td_2400"
        150   -> "td_150"
        30    -> "td_30"
    """
    sampling_rate_to_extension_map = {
        24000: "td_24k",
        2400: "td_2400",
        150: "td_150",
        30: "td_30",
    }

    try:
        return sampling_rate_to_extension_map[int(sampling_rate)]
    except (KeyError, ValueError, TypeError):
        return None


def list_unique_td_extensions(
        recording_path: str,
        subfolders=("0", "1", "2", "3", "4"),
) -> List[str]:
    """
    Scan subfolders (0..4) under a Phoenix recording folder and return unique
    file extensions that start with 'td' (e.g., 'td_24k').

    Works even if some subfolders are missing.

    Example match:
        "abc.xyz.td_24k" -> extension "td_24k"
        "data.td_24k"    -> extension "td_24k"

    """
    td_exts: Set[str] = set()

    for sf in subfolders:
        folder = os.path.join(recording_path, sf)
        if not os.path.isdir(folder):
            continue

        # listdir is faster than os.walk; switch to walk if you need recursion
        for name in os.listdir(folder):
            fp = os.path.join(folder, name)
            if not os.path.isfile(fp):
                continue

            # "extension" = text after the last dot
            base, dot, ext = name.rpartition(".")
            if dot and ext.lower().startswith("td"):
                td_exts.add(ext)  # keep original ext string (case preserved)

    return sorted(td_exts)


def get_uniform_continuous_td_counts(
        recording_path,
        file_extension: str,
        channels=("0", "1", "2", "3", "4"),
) -> Tuple[int, int]:
    counts = set()
    recording_path = pathlib.Path(recording_path)

    for ch in channels:
        ch_path = recording_path / ch

        files = sorted(ch_path.glob(f"*.{file_extension}"))
        num_files = len(files)

        if num_files == 0:
            continue

        reader = PhoenixReader.DecimatedContinuousReader(
            str(files[0]), num_files=num_files  # str() if reader expects a string path
        )

        sr = reader.header_info["sample_rate"]
        frag_s = reader.header_info["frag_period"]

        total_samples = int(sr * frag_s) * (num_files - 1)
        counts.add((num_files, total_samples))

    if not counts:
        raise ValueError("No decimated continuous files found")

    if len(counts) != 1:
        raise ValueError(f"Inconsistent TD counts across channels: {counts}")

    return counts.pop()


def get_uniform_segmented_td_counts(
        recording_path,
        file_extension: str,
        channels=("0", "1", "2", "3", "4"),
) -> Tuple[int, int]:
    """
    Returns (num_files, total_segments) for decimated *segmented* data,
    enforcing that all channels have the same counts.

    Assumes PhoenixReader.DecimatedSegmentedReader can stream across files
    using (first_file, num_files=...).
    """
    counts = set()
    recording_path = pathlib.Path(recording_path)

    for ch in channels:
        ch_path = recording_path / ch

        files = sorted(ch_path.glob(f"*.{file_extension}"))
        num_files = len(files)

        if num_files == 0:
            continue

        r = PhoenixReader.DecimatedSegmentedReader(str(files[0]), num_files=num_files)

        seg = r.read_record()

        counts.add((num_files, seg.shape[0]))

    if not counts:
        raise ValueError("No decimated segmented files found")

    if len(counts) != 1:
        raise ValueError(f"Inconsistent segmented counts across channels: {counts}")

    return counts.pop()


def read_decimated_continuous_data(
        recording_path,
        channel_map,
        file_extension
) -> Dict:
    recording_path = pathlib.Path(recording_path)
    ts = {
        'run0': {}
    }

    ts['run0']['ex'] = phoenix_readers.read_decimated_continuous(
        channel_path=recording_path / str(channel_map.get('E1')),
        file_extension=file_extension,
    )

    ts['run0']['ey'] = phoenix_readers.read_decimated_continuous(
        channel_path=recording_path / str(channel_map.get('E2')),
        file_extension=file_extension,
    )

    ts['run0']['hx'] = phoenix_readers.read_decimated_continuous(
        channel_path=recording_path / str(channel_map.get('H1')),
        file_extension=file_extension,
    )

    ts['run0']['hy'] = phoenix_readers.read_decimated_continuous(
        channel_path=recording_path / str(channel_map.get('H2')),
        file_extension=file_extension,
    )

    ts['run0']['hz'] = phoenix_readers.read_decimated_continuous(
        channel_path=recording_path / str(channel_map.get('H3')),
        file_extension=file_extension,
    )

    return ts


def read_decimated_segmented_data(recording_path, channel_map, file_extension) -> Dict:
    recording_path = pathlib.Path(recording_path)
    ts = {}

    def add_channel(chan_key: str, field: str) -> int:
        ch = channel_map.get(chan_key)
        if ch is None:
            raise KeyError(f"channel_map missing key {chan_key}")

        out_list = phoenix_readers.read_decimated_segmented(
            channel_path=recording_path / str(ch),
            file_extension=file_extension,
        )

        for num, out_val in enumerate(out_list):
            run_key = f"run{num}"
            run = ts.setdefault(run_key, {})
            run[field] = out_val["samples"]

        return len(out_list)

    counts = {
        "E1/ex": add_channel("E1", "ex"),
        "E2/ey": add_channel("E2", "ey"),
        "H1/hx": add_channel("H1", "hx"),
        "H2/hy": add_channel("H2", "hy"),
        "H3/hz": add_channel("H3", "hz"),
    }

    if len(set(counts.values())) != 1:
        raise ValueError(f"Run count mismatch across channels: {counts}")

    return ts


def optimize_time_series_dict(
        time_series_dict: dict,
        fft_length: int,
        overlap: int,
        max_runs: Optional[int] = 20,
) -> Dict:
    """
    Optimized time series dictionary
    """
    print('Optimizing time series dictionary')

    optimized_time_series_dict = {}

    for run in time_series_dict:
        optimized_time_series_dict[run] = {}
        for channel in time_series_dict[run]:
            optimized_time_series_dict[run][channel] = utils.reshape_array_with_overlap(
                window_length=fft_length,
                overlap=overlap,
                data=time_series_dict[run][channel]
            )

    del time_series_dict

    runs = list(optimized_time_series_dict.keys())
    n_runs = len(runs)

    group_size = ceil(n_runs / max_runs)

    reduced_time_series = {}

    for i in range(max_runs):
        group_runs = runs[i * group_size:(i + 1) * group_size]
        if not group_runs:
            break

        reduced_time_series[f"run{i}"] = {}

        for ch in optimized_time_series_dict[group_runs[0]].keys():
            arrays = [optimized_time_series_dict[r][ch] for r in group_runs]

            # check row consistency
            rows = {a.shape[0] for a in arrays}
            if len(rows) != 1:
                raise ValueError(f"Row mismatch in channel {ch}: {rows}")

            reduced_time_series[f"run{i}"][ch] = np.concatenate(arrays, axis=1)

    print('Optimizing time series dictionary. Done.')
    return reduced_time_series


from pathlib import Path
from collections import defaultdict


def extract_bbbbbbbb(filename: str):
    parts = filename.split("_")
    if len(parts) < 4:
        return None
    b = parts[1]
    if len(b) == 8 and all(c in "0123456789abcdefABCDEF" for c in b):
        return b.upper()
    return None


def list_files_by_bbbbbbbb_first_existing_channel(base_path, file_extension,
                                                  channels=(0, 1, 2, 3, 4)):
    """
    Search ONLY the first existing channel folder among 0..4.
    Returns: (result_dict, chosen_channel or None)
    """
    suffix = f".{file_extension}" if not file_extension.startswith(".") else file_extension
    base = Path(base_path)

    for ch in channels:
        ch_path = base / str(ch)
        if ch_path.exists() and ch_path.is_dir():
            result = defaultdict(list)
            for f in ch_path.glob(f"*{suffix}"):
                if f.is_file():
                    b = extract_bbbbbbbb(f.name)
                    if b:
                        result[b].append(f)
            return result, ch

    return defaultdict(list), None


def return_overlapping_info(file_extension, local_station_path, remote_station_path,
                            channels=(0, 1, 2, 3, 4)):
    local_map, local_ch = list_files_by_bbbbbbbb_first_existing_channel(local_station_path,
                                                                        file_extension, channels)
    remote_map, remote_ch = list_files_by_bbbbbbbb_first_existing_channel(remote_station_path,
                                                                          file_extension, channels)

    time_stamp = sorted(set(local_map.keys()) & set(remote_map.keys()))

    return time_stamp
