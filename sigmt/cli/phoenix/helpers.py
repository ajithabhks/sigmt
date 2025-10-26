import re
from pathlib import Path


def get_sampling_rates(channel_map, station_path):
    unique_sampling_rates = set()
    unique_extensions = set()

    station_path = Path(station_path)

    for ch in channel_map:
        channel_path = station_path / str(channel_map[ch])

        if not channel_path.exists():
            print(f"Missing channel folder: {channel_path}")
            continue

        # Get all file extensions in that channel folder
        extensions = {
            f.suffix.lower()
            for f in channel_path.iterdir()
            if f.is_file() and not f.name.startswith(".")
        }
        unique_extensions.update(extensions)

        # Extract sampling rates from extensions
        for e in extensions:
            match = re.search(r'td_(\d+)([kK]?)', e)
            if match:
                num = int(match.group(1))
                if match.group(2):  # if 'k' or 'K' present, multiply by 1000
                    num *= 1000
                unique_sampling_rates.add(num)

    # Sort only once after all folders are processed
    unique_sampling_rates = sorted(unique_sampling_rates)
    unique_extensions = sorted(unique_extensions)

    return unique_sampling_rates, unique_extensions
