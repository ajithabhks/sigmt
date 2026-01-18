import os
import unittest
from pathlib import Path

import numpy as np
import xarray as xr
import yaml

from sigmt.core import data_selection_tools
from sigmt.core import perform_data_selection as pds
from sigmt.core.band_averaging import BandAveraging
from sigmt.core.robust_estimation import RobustEstimation
from sigmt.utils.edi import edi_ops


def process_data(data_path: str):
    ts = {
        "run1": {
            "ex": np.load(f"{data_path}/ex.npy"),
            "ey": np.load(f"{data_path}/ey.npy"),
            "hx": np.load(f"{data_path}/hx.npy"),
            "hy": np.load(f"{data_path}/hy.npy"),
            "hz": np.load(f"{data_path}/hz.npy"),
            "rx": np.load(f"{data_path}/hx.npy"),
            "ry": np.load(f"{data_path}/hy.npy"),
        },
        "run2": {
            "ex": np.load(f"{data_path}/ex.npy"),
            "ey": np.load(f"{data_path}/ey.npy"),
            "hx": np.load(f"{data_path}/hx.npy"),
            "hy": np.load(f"{data_path}/hy.npy"),
            "hz": np.load(f"{data_path}/hz.npy"),
            "rx": np.load(f"{data_path}/hx.npy"),
            "ry": np.load(f"{data_path}/hy.npy"),
        },
    }

    with open(f"{data_path}/calibration_data_electric.yml", "r") as f:
        calibration_data_electric = yaml.safe_load(f)

    calibration_data_magnetic = {
        "instrument": "metronix",
        "hx": {
            "sensor_type": "MFS06e",
            "sensor_serial_number": 0,
            "chopper_status": 0,
            "calibration_data": {
                "chopper_off": np.load(f"{data_path}/hx_chopper_off.npy")
            },
        },
        "hy": {
            "sensor_type": "MFS06e",
            "sensor_serial_number": 0,
            "chopper_status": 0,
            "calibration_data": {
                "chopper_off": np.load(f"{data_path}/hy_chopper_off.npy")
            },
        },
        "hz": {
            "sensor_type": "MFS06e",
            "sensor_serial_number": 0,
            "chopper_status": 0,
            "calibration_data": {
                "chopper_off": np.load(f"{data_path}/hz_chopper_off.npy")
            },
        },
        "rx": {
            "sensor_type": "MFS06e",
            "sensor_serial_number": 0,
            "chopper_status": 0,
            "calibration_data": {
                "chopper_off": np.load(f"{data_path}/hx_chopper_off.npy")
            },
        },
        "ry": {
            "sensor_type": "MFS06e",
            "sensor_serial_number": 0,
            "chopper_status": 0,
            "calibration_data": {
                "chopper_off": np.load(f"{data_path}/hy_chopper_off.npy")
            },
        },
    }

    fs = 4096.0

    procinfo = {
        "localsite": "Site1",
        "processing_mode": "MT + Tipper",
        "lat": 0.0,
        "lon": 0.0,
        "start_time": 0.0,
        "elev": 0.0,
        "md_thresh": 1.5,
        "remotesite": 'Remote',
        "fs": fs,
    }

    project_setup = {
        "processing_mode": "MT + Tipper",
        "acquired_by": "COMPANY_NAME",
        "project_name": "PROJECT_NAME",
    }

    datasets = []
    for run in ts.keys():
        bandavg = BandAveraging(
            time_series=ts[run],
            sampling_frequency=procinfo["fs"],
            overlap=50,
            calibrate_electric=True,
            calibrate_magnetic=True,
            calibration_data_electric=calibration_data_electric,
            calibration_data_magnetic=calibration_data_magnetic,
            fft_length=32768,
            parzen_window_radius=0.25,
            target_frequency_table_type="default",
            frequencies_per_decade=12,
            apply_notch_filter=True,
            notch_frequency=50.0,
            process_mt=True,
            process_tipper=True,
            remote_reference=True,
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
        coh_thresh=0.8,
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
        file_name="metronix_rr",
    )


class TestMetronixRR(unittest.TestCase):
    def test_edi_matches_reference(self):
        path = Path(os.environ.get("SIGMT_TEST_DATA", r"D:/sigmt_test_data/metronix"))

        reference = path / "edi_metronix_rr_original.edi"
        self.assertTrue(reference.exists(), f"Missing reference file: {reference}")

        process_data(str(path))

        generated = path / "edi_metronix_rr.edi"
        self.assertTrue(generated.exists(), f"Generated file not found: {generated}")

        self.assertEqual(
            generated.read_bytes(),
            reference.read_bytes(),
            "Generated EDI does not match reference.",
        )


if __name__ == "__main__":
    unittest.main()
