import glob
import os
import zoneinfo
from datetime import datetime
from pathlib import Path

import numpy as np
from PhoenixGeoPy.Reader import TimeSeries as PhoenixReader


def read_decimated_continuous(
        channel_path: Path,
        file_extension: str
):
    ## DecimatedContinuousReader
    # Step 1: Get the first file (the earliest sequence)
    channel_path = Path(channel_path)
    files = sorted(glob.glob(os.path.join(channel_path, f"*.{file_extension}")))
    first_file = files[0]

    # Step 2: Tell the reader how many files to expect
    decimated_reader = PhoenixReader.DecimatedContinuousReader(first_file, num_files=len(files))

    # Step 3: Determine how many samples are in each fragment
    sample_rate = decimated_reader.header_info['sample_rate']
    frag_s = decimated_reader.header_info['frag_period']
    samples_per_file = int(sample_rate * frag_s)

    # Step 4: Read continuously until data ends
    all_data = []
    while True:
        block = decimated_reader.read_data(samples_per_file)
        if block.size == 0:  # no more data (end of last file)
            break
        all_data.append(block)

    ts = np.concatenate(all_data) * 1000  # Convert to mV

    return ts


def read_decimated_segmented(
        channel_path: Path,
        file_extension: str
):
    # Step 1: Locate all td_24K files
    channel_path = Path(channel_path)
    files = sorted(glob.glob(os.path.join(channel_path, f"*.{file_extension}")))
    first_file = files[0]

    # Step 2: Create a reader that can stream across them
    r = PhoenixReader.DecimatedSegmentedReader(first_file, num_files=len(files))

    # Step 3: Read every segment across all files
    segments = []

    while True:
        seg = r.read_record()  # read next segment (data + subheader)
        if seg.size == 0:  # end of all files
            break
        segments.append({
            "timestamp": r.subheader["timestamp"] * 1000,  # to mV
            "samples": seg,
            "minVal": r.subheader["minVal"],
            "maxVal": r.subheader["maxVal"],
            "avgVal": r.subheader["avgVal"],
        })

    print(f"Read {len(segments)} segments across {len(files)} files")

    return segments
