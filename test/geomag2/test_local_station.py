import os
import unittest
from pathlib import Path

import numpy as np
import xarray as xr
import yaml
from scipy import signal

from sigmt.core import data_selection_tools
from sigmt.core import perform_data_selection as pds
from sigmt.core import plots
from sigmt.core.band_averaging import BandAveraging
from sigmt.core.robust_estimation import RobustEstimation
from sigmt.utils.edi import edi_ops
from sigmt.utils.geomag import read_time_series as geomag_reader


def process_data(data_path: str):
    header, df = geomag_reader.read_geomag_ascii(rf"{data_path}\MS_12_190715060407.TXT")

    time_series = {
        "run1": {
            "ex": df['ex'].values,
            "ey": df['ey'].values,
            "hx": df['hx'].values,
            "hy": df['hy'].values,
            "hz": df['hz'].values,
        }
    }

    fs = 1.0

    decimation_factors = [4, 4]
    for decimation_factor in decimation_factors:
        for run in time_series:
            for channel in time_series[run]:
                time_series[run][channel] = signal.decimate(
                    time_series[run][channel],
                    decimation_factor,
                    n=None,
                    ftype='iir'
                )
        fs = fs / decimation_factor

    with open(f"{data_path}/calibration_data_electric.yml", "r") as f:
        calibration_data_electric = yaml.safe_load(f)

    procinfo = {
        "localsite": "Site1",
        "processing_mode": "MT + Tipper",
        "lat": 0.0,
        "lon": 0.0,
        "start_time": 0.0,
        "elev": 0.0,
        "md_thresh": 1.5,
        "remotesite": None,
        "fs": fs,
    }

    project_setup = {
        "processing_mode": "MT + Tipper",
        "acquired_by": "COMPANY_NAME",
        "project_name": "PROJECT_NAME",
    }

    datasets = []
    for run in time_series.keys():
        bandavg = BandAveraging(
            time_series=time_series[run],
            sampling_frequency=procinfo["fs"],
            overlap=50,
            calibrate_electric=True,
            calibrate_magnetic=False,
            calibration_data_electric=calibration_data_electric,
            fft_length=4096,
            parzen_window_radius=0.6,
            target_frequency_table_type="default",
            frequencies_per_decade=12,
            process_mt=True,
            process_tipper=True,
            remote_reference=False,
        )
        datasets.append(bandavg.band_averaged_dataset)

    bandavg_dataset = xr.concat(datasets, dim="time_window").assign_coords(
        time_window=np.arange(len(xr.concat(datasets, dim="time_window").time_window))
    )

    bandavg_dataset["dof"] = xr.DataArray(
        bandavg.dof,
        coords={"frequency": bandavg_dataset.coords["frequency"]},
        dims="frequency",
    )

    bandavg_dataset["coh_ex"] = (
        ("time_window", "frequency"),
        data_selection_tools.cohex(bandavg_dataset),
    )
    bandavg_dataset["coh_ey"] = (
        ("time_window", "frequency"),
        data_selection_tools.cohey(bandavg_dataset),
    )
    bandavg_dataset["coh_hz"] = (
        ("time_window", "frequency"),
        data_selection_tools.cohhz(bandavg_dataset),
    )
    bandavg_dataset["alpha_h"], bandavg_dataset["alpha_e"] = data_selection_tools.pdvalues(
        bandavg_dataset
    )

    bandavg_dataset = pds.perform_coh_thresh(
        bandavg_dataset,
        coh_thresh=0.95,
        min_percent=20,
        channels=["ex", "ey", "hz"],
    )

    estimation_instance = RobustEstimation(procinfo, bandavg_dataset)
    estimates = estimation_instance.estimates

    edi_ops.save_edi(
        estimates=estimates,
        procinfo=procinfo,
        project_setup=project_setup,
        save_path=data_path,
        file_name="metronix_local",
    )

    plots.plot_mt_app_res(
        dataset=estimates,
        procinfo=procinfo,
        save_path=data_path
    )

    plots.plot_tipper(
        dataset=estimates,
        procinfo=procinfo,
        save_path=data_path
    )

    plots.plot_coherency(
        dataset=estimates,
        procinfo=procinfo,
        save_path=data_path
    )


class TestMetronixLocal(unittest.TestCase):
    def test_geomag_processing(self):
        path = Path(os.environ.get("SIGMT_TEST_DATA", r"D:/MTDATA/geomag2"))

        process_data(str(path))


if __name__ == "__main__":
    unittest.main()
