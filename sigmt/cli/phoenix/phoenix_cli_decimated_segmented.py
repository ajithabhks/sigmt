import json
from pathlib import Path

import numpy as np
import xarray as xr

import sigmt.cli.phoenix.data_readers as phoenix_readers
from sigmt.core import data_selection_tools, plots
from sigmt.core import perform_data_selection as pds
from sigmt.core.band_averaging import BandAveraging
from sigmt.core.robust_estimation import RobustEstimation
from sigmt.utils.edi import edi_ops

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

fs = 24000
fft_length = 16384
parzen_window_radius = 0.25
notch_frequency = 60
overlap = 50

procinfo = {}
procinfo["localsite"] = "MT1"
procinfo['processing_mode'] = "MT + Tipper"
procinfo['lat'] = recmeta_data['timing']['gps_lat']
procinfo['lon'] = recmeta_data['timing']['gps_lon']
procinfo['elev'] = recmeta_data['timing']['gps_alt']
procinfo['start_time'] = 0

procinfo['md_thresh'] = 1.5

procinfo['remotesite'] = None
procinfo["fs"] = fs

project_setup = {}
project_setup['processing_mode'] = 'MT + Tipper'
project_setup['acquired_by'] = 'COMPANY_NAME'
project_setup['project_name'] = 'PROJECT_NAME'

save_path = "D:/"

hx_coil_serial = recmeta_data['chconfig']['chans'][channel_map['H1']]['serial']
hy_coil_serial = recmeta_data['chconfig']['chans'][channel_map['H2']]['serial']
hz_coil_serial = recmeta_data['chconfig']['chans'][channel_map['H3']]['serial']

cal_data_path = Path('D:/MTDATA/PhoenixData/Cals')
hx_coil_data_path = cal_data_path / f'{hx_coil_serial}.json'
with open(hx_coil_data_path, 'r', encoding='utf-8') as f:
    hx_coil_data = json.load(f)
hy_coil_data_path = cal_data_path / f'{hy_coil_serial}.json'
with open(hx_coil_data_path, 'r', encoding='utf-8') as f:
    hy_coil_data = json.load(f)
hz_coil_data_path = cal_data_path / f'{hz_coil_serial}.json'
with open(hx_coil_data_path, 'r', encoding='utf-8') as f:
    hz_coil_data = json.load(f)

calibration_data_magnetic = {
    'instrument': 'phoenix',
    'hx': {},
    'hy': {},
    'hz': {}
}

calibration_data_electric = {
    'ex': {},
    'ey': {},
}

calibration_data_electric['ex']['x1'] = abs(recmeta_data['chconfig']['chans']
                                            [channel_map['E1']]['length1'])
calibration_data_electric['ex']['x2'] = abs(recmeta_data['chconfig']['chans']
                                            [channel_map['E1']]['length2'])
calibration_data_electric['ey']['y1'] = abs(recmeta_data['chconfig']['chans']
                                            [channel_map['E2']]['length1'])
calibration_data_electric['ey']['y2'] = abs(recmeta_data['chconfig']['chans']
                                            [channel_map['E2']]['length2'])

calibration_data_magnetic['hx']['calibration_data'] = hx_coil_data['cal_data'][0]['chan_data'][0]
calibration_data_magnetic['hy']['calibration_data'] = hy_coil_data['cal_data'][0]['chan_data'][0]
calibration_data_magnetic['hz']['calibration_data'] = hz_coil_data['cal_data'][0]['chan_data'][0]

ts = {}

out_list = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('E1')),
    file_extension="td_24k",
)
for num, out_val in enumerate(out_list):
    ts[f'run{num}'] = {}
    ts[f'run{num}']['ex'] = out_val['samples']

out_list = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('E2')),
    file_extension="td_24k",
)
for num, out_val in enumerate(out_list):
    ts[f'run{num}']['ey'] = out_val['samples']

out_list = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('H1')),
    file_extension="td_24k",
)
for num, out_val in enumerate(out_list):
    ts[f'run{num}']['hx'] = out_val['samples']

out_list = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('H2')),
    file_extension="td_24k",
)
for num, out_val in enumerate(out_list):
    ts[f'run{num}']['hy'] = out_val['samples']

out_list = phoenix_readers.read_decimated_segmented(
    channel_path=station_path / str(channel_map.get('H3')),
    file_extension="td_24k",
)
for num, out_val in enumerate(out_list):
    ts[f'run{num}']['hz'] = out_val['samples']

# import scipy.signal as signal
# ts['run1']['ex'] = signal.decimate(ts['run1']['ex'], 4, n=None, ftype='iir')
# ts['run1']['ey'] = signal.decimate(ts['run1']['ey'], 4, n=None, ftype='iir')
# ts['run1']['hx'] = signal.decimate(ts['run1']['hx'], 4, n=None, ftype='iir')
# ts['run1']['hy'] = signal.decimate(ts['run1']['hy'], 4, n=None, ftype='iir')
# ts['run1']['hz'] = signal.decimate(ts['run1']['hz'], 4, n=None, ftype='iir')
# procinfo["fs"] = procinfo["fs"]/4

datasets = []
for run in ts.keys():
    bandavg = BandAveraging(
        time_series=ts[run],
        sampling_frequency=fs,
        fft_length=fft_length,
        parzen_window_radius=parzen_window_radius,
        overlap=overlap,
        frequencies_per_decade=12,
        remote_reference=False,
        calibrate_electric=True,
        calibrate_magnetic=True,
        calibration_data_electric=calibration_data_electric,
        calibration_data_magnetic=calibration_data_magnetic,
        instrument='phoenix',
        apply_notch_filter=False,
        notch_frequency=notch_frequency,
        process_mt=True,
        process_tipper=True,
    )
    datasets.append(bandavg.band_averaged_dataset)

bandavg_dataset = xr.concat(
    datasets, dim='time_window'
).assign_coords(
    time_window=np.arange(
        len(xr.concat(
            datasets,
            dim='time_window'
        ).time_window
            )
    )
)
dof = bandavg.dof
avgf = bandavg.avgf
bandavg_dataset['dof'] = xr.DataArray(
    dof,
    coords={'frequency': bandavg_dataset.coords['frequency']},
    dims='frequency'
)

bandavg_dataset['coh_ex'] = (
    ('time_window', 'frequency'), data_selection_tools.cohex(bandavg_dataset))
bandavg_dataset['coh_ey'] = (
    ('time_window', 'frequency'), data_selection_tools.cohey(bandavg_dataset))
bandavg_dataset['coh_hz'] = (
    ('time_window', 'frequency'), data_selection_tools.cohhz(bandavg_dataset))
bandavg_dataset['alpha_h'], bandavg_dataset['alpha_e'] = data_selection_tools.pdvalues(
    bandavg_dataset)

# Apply coherency threshold
bandavg_dataset = pds.perform_coh_thresh(
    bandavg_dataset,
    coh_thresh=0.8,
    min_percent=20,
    channels=['ex', 'ey', 'hz']
)

estimation_instance = RobustEstimation(procinfo, bandavg_dataset)
estimates = estimation_instance.estimates

plots.plot_mt_app_res(
    dataset=estimates,
    procinfo=procinfo,
    save_path=save_path
)
plots.plot_tipper(
    dataset=estimates,
    procinfo=procinfo,
    save_path=save_path
)
plots.plot_coherency(
    dataset=estimates,
    procinfo=procinfo,
    save_path=save_path
)
edi_ops.save_edi(
    estimates=estimates,
    procinfo=procinfo,
    project_setup=project_setup,
    save_path=save_path
)

"""
Example Unique
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
"""
