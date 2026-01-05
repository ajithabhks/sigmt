"""
Utility functions for the phoenix
"""

import os
import pathlib
import re
from typing import List, Set, Tuple, Dict

from PhoenixGeoPy.Reader import TimeSeries as PhoenixReader

import sigmt.cli.phoenix.data_readers as phoenix_readers


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

        total_samples = int(sr * frag_s) * (num_files-1)
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


    # temp = np.hstack(temp)
    # total_cols = temp.shape[1]
    #
    # for i, start in enumerate(range(0, total_cols, max_cols), start=1):
    #     end = start + max_cols
    #     ts.setdefault(f'run{i}', {})['ex'] = temp[:, start:end]
    #
    # out_list = phoenix_readers.read_decimated_segmented(
    #     channel_path=station_path / str(channel_map.get('E2')),
    #     file_extension="td_24k",
    # )
    # temp = []
    # for num, out_val in enumerate(out_list):
    #     temp.append(utils.reshape_array_with_overlap(
    #         window_length=fft_length,
    #         overlap=overlap,
    #         data=out_val['samples']
    #     ))
    # temp = np.hstack(temp)
    # total_cols = temp.shape[1]
    #
    # for i, start in enumerate(range(0, total_cols, max_cols), start=1):
    #     end = start + max_cols
    #     ts.setdefault(f'run{i}', {})['ey'] = temp[:, start:end]
    #
    # out_list = phoenix_readers.read_decimated_segmented(
    #     channel_path=station_path / str(channel_map.get('H1')),
    #     file_extension="td_24k",
    # )
    # temp = []
    # for num, out_val in enumerate(out_list):
    #     temp.append(utils.reshape_array_with_overlap(
    #         window_length=fft_length,
    #         overlap=overlap,
    #         data=out_val['samples']
    #     ))
    # temp = np.hstack(temp)
    # total_cols = temp.shape[1]
    #
    # for i, start in enumerate(range(0, total_cols, max_cols), start=1):
    #     end = start + max_cols
    #     ts.setdefault(f'run{i}', {})['hx'] = temp[:, start:end]
    #
    # out_list = phoenix_readers.read_decimated_segmented(
    #     channel_path=station_path / str(channel_map.get('H2')),
    #     file_extension="td_24k",
    # )
    # temp = []
    # for num, out_val in enumerate(out_list):
    #     temp.append(utils.reshape_array_with_overlap(
    #         window_length=fft_length,
    #         overlap=overlap,
    #         data=out_val['samples']
    #     ))
    # temp = np.hstack(temp)
    # total_cols = temp.shape[1]
    #
    # for i, start in enumerate(range(0, total_cols, max_cols), start=1):
    #     end = start + max_cols
    #     ts.setdefault(f'run{i}', {})['hy'] = temp[:, start:end]
    #
    # out_list = phoenix_readers.read_decimated_segmented(
    #     channel_path=station_path / str(channel_map.get('H3')),
    #     file_extension="td_24k",
    # )
    # temp = []
    # for num, out_val in enumerate(out_list):
    #     temp.append(utils.reshape_array_with_overlap(
    #         window_length=fft_length,
    #         overlap=overlap,
    #         data=out_val['samples']
    #     ))
    # temp = np.hstack(temp)
    # total_cols = temp.shape[1]
    #
    # for i, start in enumerate(range(0, total_cols, max_cols), start=1):
    #     end = start + max_cols
    #     ts.setdefault(f'run{i}', {})['hz'] = temp[:, start:end]
    #
    # return ts