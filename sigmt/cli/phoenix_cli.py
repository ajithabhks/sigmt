from pathlib import Path
import json
from sigmt.utils.phoenix.phoenix_reader import PhoenixMTReader

station_path = Path(r'D:\MTDATA\PhoenixData\SampleData1\10771_2024-08-26-152914')

# Define the file path
recmeta_file = station_path / 'recmeta.json'

# Read JSON file
with open(recmeta_file, 'r', encoding='utf-8') as f:
    recmeta_data = json.load(f)

channel_map = recmeta_data['channel_map']['mapping']
procinfo = {}
procinfo['lat'] = recmeta_data['timing']['gps_lat']
procinfo['lon'] = recmeta_data['timing']['gps_lon']
procinfo['elev'] = recmeta_data['timing']['gps_alt']


header = {}

recording_path = station_path / '0' / '10771_66CC9F4A_0_0000000A.td_150'
# Works only for 150 s/s
reader = PhoenixMTReader(recording_path)
header, data = reader.read()
pass
