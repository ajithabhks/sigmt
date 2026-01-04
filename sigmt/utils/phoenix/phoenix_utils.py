"""
Utility functions for the phoenix
"""

import os
import re
from typing import List, Set


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
