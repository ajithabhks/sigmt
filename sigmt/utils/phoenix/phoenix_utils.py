"""
Utility functions for the phoenix
"""

import json
import os
import pathlib
import re
from collections import defaultdict
from datetime import datetime, timedelta
from math import ceil
from pathlib import Path
from typing import List, Set, Dict, Optional, Any, Tuple

import numpy as np

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
    :type project_dir: str
    :return: List of valid site names
    :rtype: List
    """

    # allow anything after the timestamp
    recording_pattern = re.compile(
        r"^\d+_\d{4}-\d{2}-\d{2}-\d{6}(?:_.+)?$"
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


def get_sampling_rate_list(recording_path: str) -> List[str]:
    """
    Return a sorted list of unique sampling rates found in a recording folder.

    Mapping:
        td_24k  -> 24000
        td_2400 -> 2400
        td_150  -> 150
        td_30   -> 30

    :param recording_path: Path to the recording folder.
    :type recording_path: str
    :return: List of sampling rates
    :rtype: List
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


def sampling_rate_to_extension(sampling_rate: int) -> str:
    """
    Convert a sampling rate to its corresponding td_* extension.

    Examples:
        24000 -> "td_24k"
        2400  -> "td_2400"
        150   -> "td_150"
        30    -> "td_30"

    :param sampling_rate: Sampling rate
    :type sampling_rate: int
    :return: File extension
    :rtype: str
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
        subfolders: tuple = ("0", "1", "2", "3", "4"),
) -> List[str]:
    """
    Scan subfolders (0..4) under a Phoenix recording folder and return unique
    file extensions that start with 'td' (e.g., 'td_24k').

    Works even if some subfolders are missing.

    Example match:
        "abc.xyz.td_24k" -> extension "td_24k"
        "data.td_24k"    -> extension "td_24k"

    :param recording_path: Path to recording path
    :type recording_path: str
    :param subfolders: Folders in recording path
    :type subfolders: tuple
    :return: List of file extensions
    :rtype: List
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


def read_decimated_continuous_data(
        recording_path: str,
        channel_map: Dict,
        file_extension: str
) -> Tuple[Dict, Optional[datetime]]:
    """
    Read decimated continuous data.

    :param recording_path: Path to the recording folder
    :type recording_path: str
    :param channel_map: Channel map dict
    :type channel_map: Dict
    :param file_extension: File extension of file
    :type file_extension: str
    :return: Dictionary containing time series data and timestamp
    :rtype: Tuple

    """
    recording_path = pathlib.Path(recording_path)

    ts = {
        "run0": {}
    }

    timestamps = {}

    channel_keys = {
        "E1": "ex",
        "E2": "ey",
        "H1": "hx",
        "H2": "hy",
        "H3": "hz",
    }

    for map_key, component in channel_keys.items():
        channel = channel_map.get(map_key)

        if channel is None:
            continue

        channel_path = recording_path / str(channel)

        if not channel_path.exists():
            print(f"Warning: channel path does not exist: {channel_path}")
            continue

        data, timestamp = phoenix_readers.read_decimated_continuous(
            channel_path=channel_path,
            file_extension=file_extension,
        )

        ts["run0"][component] = data
        timestamps[component] = timestamp

    if timestamps:
        unique_timestamps = set(timestamps.values())

        if len(unique_timestamps) != 1:
            raise ValueError(
                f"Channel timestamps do not match: {timestamps}"
            )
        timestamp = next(iter(unique_timestamps))
    else:
        timestamp = None

    return ts, timestamp


def read_decimated_segmented_data(
        recording_path: str,
        channel_map: Dict,
        file_extension: str,
) -> Dict:
    """
    Read decimated segmented data

    :param recording_path: Path to recordings folder
    :type recording_path: str
    :param channel_map: Dict of channel map
    :type channel_map: dict
    :param file_extension: File extension
    :type file_extension: str
    :return: Dict of time series
    :rtype: dict

    """
    recording_path = pathlib.Path(recording_path)
    ts: Dict[str, Dict[str, Any]] = {}

    ex_map = channel_map.get('E1', None)
    if ex_map is not None:
        channel_path = recording_path / str(ex_map)
        if channel_path.exists():
            segments = phoenix_readers.read_decimated_segmented(
                channel_path=channel_path,
                file_extension=file_extension,
            )
            for num, list_item in enumerate(segments):
                run_key = f"run{num}"
                current_timestamp = list_item.get("timestamp")
                run_data = ts.setdefault(run_key, {})
                # Check timestamp if it was already set
                if "timestamp" in run_data:
                    if run_data["timestamp"] != current_timestamp:
                        raise ValueError(
                            f"Timestamp mismatch for {run_key}: "
                            f"existing={run_data['timestamp']}, "
                            f"current={current_timestamp}"
                        )
                else:
                    run_data["timestamp"] = current_timestamp
                run_data["ex"] = list_item.get("samples")

    ey_map = channel_map.get('E2', None)
    if ey_map is not None:
        channel_path = recording_path / str(ey_map)
        if channel_path.exists():
            segments = phoenix_readers.read_decimated_segmented(
                channel_path=channel_path,
                file_extension=file_extension,
            )
            for num, list_item in enumerate(segments):
                run_key = f"run{num}"
                current_timestamp = list_item.get("timestamp")
                run_data = ts.setdefault(run_key, {})
                # Check timestamp if it was already set
                if "timestamp" in run_data:
                    if run_data["timestamp"] != current_timestamp:
                        raise ValueError(
                            f"Timestamp mismatch for {run_key}: "
                            f"existing={run_data['timestamp']}, "
                            f"current={current_timestamp}"
                        )
                else:
                    run_data["timestamp"] = current_timestamp
                run_data["ey"] = list_item.get("samples")

    hx_map = channel_map.get('H1', None)
    if hx_map is not None:
        channel_path = recording_path / str(hx_map)
        if channel_path.exists():
            segments = phoenix_readers.read_decimated_segmented(
                channel_path=channel_path,
                file_extension=file_extension,
            )
            for num, list_item in enumerate(segments):
                run_key = f"run{num}"
                current_timestamp = list_item.get("timestamp")
                run_data = ts.setdefault(run_key, {})
                # Check timestamp if it was already set
                if "timestamp" in run_data:
                    if run_data["timestamp"] != current_timestamp:
                        raise ValueError(
                            f"Timestamp mismatch for {run_key}: "
                            f"existing={run_data['timestamp']}, "
                            f"current={current_timestamp}"
                        )
                else:
                    run_data["timestamp"] = current_timestamp
                run_data["hx"] = list_item.get("samples")

    hy_map = channel_map.get('H2', None)
    if hy_map is not None:
        channel_path = recording_path / str(hy_map)
        if channel_path.exists():
            segments = phoenix_readers.read_decimated_segmented(
                channel_path=channel_path,
                file_extension=file_extension,
            )
            for num, list_item in enumerate(segments):
                run_key = f"run{num}"
                current_timestamp = list_item.get("timestamp")
                run_data = ts.setdefault(run_key, {})
                # Check timestamp if it was already set
                if "timestamp" in run_data:
                    if run_data["timestamp"] != current_timestamp:
                        raise ValueError(
                            f"Timestamp mismatch for {run_key}: "
                            f"existing={run_data['timestamp']}, "
                            f"current={current_timestamp}"
                        )
                else:
                    run_data["timestamp"] = current_timestamp
                run_data["hy"] = list_item.get("samples")

    hz_map = channel_map.get('H3', None)
    if hz_map is not None:
        channel_path = recording_path / str(hz_map)
        if channel_path.exists():
            segments = phoenix_readers.read_decimated_segmented(
                channel_path=channel_path,
                file_extension=file_extension,
            )
            for num, list_item in enumerate(segments):
                run_key = f"run{num}"
                current_timestamp = list_item.get("timestamp")
                run_data = ts.setdefault(run_key, {})
                # Check timestamp if it was already set
                if "timestamp" in run_data:
                    if run_data["timestamp"] != current_timestamp:
                        raise ValueError(
                            f"Timestamp mismatch for {run_key}: "
                            f"existing={run_data['timestamp']}, "
                            f"current={current_timestamp}"
                        )
                else:
                    run_data["timestamp"] = current_timestamp
                run_data["hz"] = list_item.get("samples")
    return ts


def optimize_time_series_dict(
        time_series_dict: dict,
        fft_length: int,
        overlap: int,
        max_runs: Optional[int] = 20,
) -> Dict:
    """
    Optimized time series dictionary. It reduces number of runs.

    :param time_series_dict: Dictionary of time series
    :type time_series_dict: dict
    :param fft_length: FFT Length
    :type fft_length: int
    :param overlap: Overlap in percentage
    :type overlap: int
    :param max_runs: Number of runs reduced to
    :type max_runs: int
    :return: Optimized time series dict
    :rtype: dict
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


def extract_bbbbbbbb(filename: str) -> str:
    """
    Extract time stamp from filename
    :param filename: filename
    :type filename: str
    :return: string
    :rtype: str
    """
    parts = filename.split("_")
    if len(parts) < 4:
        return None
    b = parts[1]
    if len(b) == 8 and all(c in "0123456789abcdefABCDEF" for c in b):
        return b.upper()
    return None


def list_files_by_bbbbbbbb_first_existing_channel(
        base_path: str,
        file_extension: str,
        channels: tuple = (0, 1, 2, 3, 4)
) -> tuple:
    """
    Search ONLY the first existing channel folder among 0..4.
    Returns: (result_dict, chosen_channel or None)

    :param base_path: Path to folder
    :type base_path: str
    :param file_extension: File extension
    :type file_extension: str
    :param channels: Channels
    :type channels: tuple

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


def return_overlapping_info(
        file_extension,
        local_station_path,
        remote_station_path,
        channels=(0, 1, 2, 3, 4)
):
    local_map, local_ch = list_files_by_bbbbbbbb_first_existing_channel(local_station_path,
                                                                        file_extension, channels)
    remote_map, remote_ch = list_files_by_bbbbbbbb_first_existing_channel(remote_station_path,
                                                                          file_extension, channels)

    time_stamp = sorted(set(local_map.keys()) & set(remote_map.keys()))

    return time_stamp


def prepare_calibration_data_electric(
        local_recmeta_data: Dict,
        channel_map: Dict
) -> Dict:
    """
    Prepare electric calibration data.
    """
    calibration_data_electric = {
        'ex': {},
        'ey': {},
    }
    calibration_data_electric['ex']['x1'] = np.abs(
        local_recmeta_data['chconfig']['chans']
        [channel_map['E1']]['length1']
    )
    calibration_data_electric['ex']['x2'] = np.abs(
        local_recmeta_data['chconfig']['chans']
        [channel_map['E1']]['length2']
    )
    calibration_data_electric['ey']['y1'] = np.abs(
        local_recmeta_data['chconfig']['chans']
        [channel_map['E2']]['length1']
    )
    calibration_data_electric['ey']['y2'] = np.abs(
        local_recmeta_data['chconfig']['chans']
        [channel_map['E2']]['length2']
    )

    return calibration_data_electric


def prepare_calibration_data_magnetic(
        project_dir,
        local_recmeta_data: Dict[str, Any],
        local_channel_map: Dict[str, Any],
        remote_recmeta_data: Optional[Dict[str, Any]] = None,
        remote_channel_map: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Prepare magnetic calibration data for Phoenix instrument.

    Expects calibration files at: <project_dir>/calibration_files/<serial>.json
    Extracts: data['cal_data'][0]['chan_data'][0]
    """

    cal_dir = pathlib.Path(project_dir) / "calibration_files"

    def _extract_chan_cal(serial: str) -> Any:
        path = cal_dir / f"{serial}.json"
        try:
            with path.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Missing calibration file: {path}") from e

        # Keep same extraction semantics, but with a clearer error if shape differs
        try:
            return data["cal_data"][0]["chan_data"][0]
        except (KeyError, IndexError, TypeError) as e:
            raise ValueError(
                f"Unexpected calibration JSON structure in {path}. "
                f"Expected data['cal_data'][0]['chan_data'][0]."
            ) from e

    def _serial_from(meta: Dict[str, Any], chan_map: Dict[str, Any], key: str) -> str:
        try:
            chan_idx = chan_map[key]
            return meta["chconfig"]["chans"][chan_idx]["serial"]
        except KeyError as e:
            raise KeyError(f"Missing key while resolving serial for {key}: {e}") from e
        except (IndexError, TypeError) as e:
            raise ValueError(f"Bad channel map/index while resolving serial for {key}.") from e

    calibration_data_magnetic: Dict[str, Any] = {
        "instrument": "phoenix",
        "hx": {"calibration_data": _extract_chan_cal(
            _serial_from(local_recmeta_data, local_channel_map, "H1"))},
        "hy": {"calibration_data": _extract_chan_cal(
            _serial_from(local_recmeta_data, local_channel_map, "H2"))},
        "hz": {"calibration_data": _extract_chan_cal(
            _serial_from(local_recmeta_data, local_channel_map, "H3"))},
        "rx": {},
        "ry": {},
    }

    if remote_recmeta_data and remote_channel_map:
        calibration_data_magnetic["rx"]["calibration_data"] = _extract_chan_cal(
            _serial_from(remote_recmeta_data, remote_channel_map, "H1")
        )
        calibration_data_magnetic["ry"]["calibration_data"] = _extract_chan_cal(
            _serial_from(remote_recmeta_data, remote_channel_map, "H2")
        )

    return calibration_data_magnetic


def get_time_overlap(local_start, local_stop, remote_start, remote_stop):
    """
    Find overlap between two Unix timestamp ranges.

    Returns:
        (overlap_seconds, overlap_hms)
        or
        (0, "00:00:00") if there is no overlap.
    """
    if None in (local_start, local_stop, remote_start, remote_stop):
        return 0, "00:00:00"

    overlap_start = max(local_start, remote_start)
    overlap_stop = min(local_stop, remote_stop)

    overlap_seconds = max(0, overlap_stop - overlap_start)

    hours, remainder = divmod(overlap_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    overlap_hms = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

    return overlap_seconds, overlap_hms


def count_matching_timestamps(local_ts, remote_ts):
    local_timestamps = {
        data.get("timestamp")
        for data in local_ts.values()
        if data.get("timestamp") is not None
    }

    remote_timestamps = {
        data.get("timestamp")
        for data in remote_ts.values()
        if data.get("timestamp") is not None
    }

    matching_timestamps = local_timestamps & remote_timestamps

    return len(matching_timestamps), matching_timestamps


def build_timestamp_map(time_series, name):
    """Create timestamp -> run data mapping."""

    timestamp_map = {}

    for run_key, run_data in time_series.items():
        timestamp = run_data.get("timestamp")

        if timestamp is None:
            raise ValueError(
                f"{name}: timestamp missing in {run_key}."
            )

        if timestamp in timestamp_map:
            raise ValueError(
                f"{name}: duplicate timestamp {timestamp} found."
            )

        timestamp_map[timestamp] = run_data

    return timestamp_map


def get_array_length(run_data, name, timestamp):
    """Check that all arrays in a run have the same length."""

    lengths = {}

    for channel, data in run_data.items():

        if channel == "timestamp":
            continue

        if isinstance(data, np.ndarray):
            if data.ndim != 1:
                raise ValueError(
                    f"{name}: {channel} at timestamp {timestamp} "
                    f"has invalid shape {data.shape}."
                )

            lengths[channel] = len(data)

    if not lengths:
        raise ValueError(
            f"{name}: no channel arrays found at timestamp {timestamp}."
        )

    if len(set(lengths.values())) != 1:
        raise ValueError(
            f"{name}: inconsistent array lengths at timestamp "
            f"{timestamp}: {lengths}"
        )

    return next(iter(lengths.values()))


def trim_to_matching_timestamps(local_ts, remote_ts):
    """
    Keep only matching timestamps, sort in ascending order,
    verify local/remote sample counts, and remove timestamp
    from the final dictionaries.
    """

    local_map = build_timestamp_map(local_ts, "Local")
    remote_map = build_timestamp_map(remote_ts, "Remote")

    matching_timestamps = sorted(
        set(local_map) & set(remote_map)
    )

    if not matching_timestamps:
        raise ValueError(
            "No matching timestamps found between local and remote data."
        )

    local_trimmed = {}
    remote_trimmed = {}

    for num, timestamp in enumerate(matching_timestamps):

        local_data = local_map[timestamp]
        remote_data = remote_map[timestamp]

        local_length = get_array_length(
            local_data,
            "Local",
            timestamp
        )

        remote_length = get_array_length(
            remote_data,
            "Remote",
            timestamp
        )

        if local_length != remote_length:
            raise ValueError(
                f"Sample count mismatch at timestamp {timestamp}: "
                f"local={local_length}, remote={remote_length}."
            )

        run_key = f"run{num}"

        local_trimmed[run_key] = {
            key: value
            for key, value in local_data.items()
            if key != "timestamp"
        }

        remote_trimmed[run_key] = {
            key: value
            for key, value in remote_data.items()
            if key != "timestamp"
        }

    return local_trimmed, remote_trimmed


def align_continuous_time_series(
        local_ts: Dict,
        local_timestamp: datetime,
        remote_ts: Dict,
        remote_timestamp: datetime,
        sampling_rate: float,
) -> Tuple[Dict, Dict]:
    """
    Trim local and remote time series to their overlapping time interval.

    Assumes all channels within each time series have the same number of samples.
    Timestamps are datetime objects.
    """
    if sampling_rate <= 0:
        raise ValueError("Sampling rate must be greater than zero.")

    if not local_ts.get("run0") or not remote_ts.get("run0"):
        raise ValueError("Local or remote time series is empty.")

    local_n = min(len(data) for data in local_ts["run0"].values())
    remote_n = min(len(data) for data in remote_ts["run0"].values())

    local_end = local_timestamp + timedelta(
        seconds=local_n / sampling_rate
    )
    remote_end = remote_timestamp + timedelta(
        seconds=remote_n / sampling_rate
    )

    overlap_start = max(local_timestamp, remote_timestamp)
    overlap_end = min(local_end, remote_end)

    if overlap_start >= overlap_end:
        raise ValueError("Local and remote time series do not overlap.")

    local_start_idx = round(
        (overlap_start - local_timestamp).total_seconds()
        * sampling_rate
    )

    remote_start_idx = round(
        (overlap_start - remote_timestamp).total_seconds()
        * sampling_rate
    )

    overlap_n = round(
        (overlap_end - overlap_start).total_seconds()
        * sampling_rate
    )

    if overlap_n <= 0:
        raise ValueError("No valid samples exist within the overlap.")

    for channel in local_ts["run0"]:
        local_ts["run0"][channel] = local_ts["run0"][channel][
            local_start_idx: local_start_idx + overlap_n
        ]

    for channel in remote_ts["run0"]:
        remote_ts["run0"][channel] = remote_ts["run0"][channel][
            remote_start_idx: remote_start_idx + overlap_n
        ]

    return local_ts, remote_ts
