import json
from pathlib import Path

import sigmt.cli.phoenix.data_readers as phoenix_readers
import sigmt.cli.phoenix.helpers as phoenix_helpers

station_path = Path(r'D:\MTDATA\PhoenixData\SampleData1\10771_2024-08-26-152914')

# Define the file path
recmeta_file = station_path / 'recmeta.json'

# Read JSON file
with open(recmeta_file, 'r', encoding='utf-8') as f:
    recmeta_data = json.load(f)

channel_map = recmeta_data['channel_map']['mapping']
channel_map = {ch['tag']: ch['idx'] for ch in channel_map}
ex = channel_map.get('E1')
ey = channel_map.get('E2')
hx = channel_map.get('H1')
hy = channel_map.get('H2')
hz = channel_map.get('H3')

procinfo = {}
procinfo['lat'] = recmeta_data['timing']['gps_lat']
procinfo['lon'] = recmeta_data['timing']['gps_lon']
procinfo['elev'] = recmeta_data['timing']['gps_alt']

folder = station_path / str(channel_map.get('E1'))
sampling_rates, file_extensions = phoenix_helpers.get_sampling_rates(channel_map, station_path)

# DecimatedContinuous
ex_data1 = phoenix_readers.read_decimated_continuous(
    channel_path=station_path / str(channel_map.get('E1')),
    file_extension="td_150",
)

# DecimatedSegmented
ex_data2 = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('E1')),
    file_extension="td_24K",
)
