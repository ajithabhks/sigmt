import numpy as np


def read_calibration_metronix_txt(filepath):
    """
    Read calibration data "chopper_on" and "chopper_off" from Metronix calibration
    txt file.

    :param filepath: File path
    :type filepath: str or Path

    :return: Dictionary containing chopper_off and chopper_on data as numpy arrays
    :rtype: dict

    """
    calibration_data = {}
    with open(filepath, "r") as f:
        lines = f.readlines()

    chopper_on = []
    chopper_off = []

    current_section = None

    for line in lines:
        line = line.strip()

        # Detect section headers
        if line.startswith("Chopper On"):
            current_section = "on"
            continue
        elif line.startswith("Chopper Off"):
            current_section = "off"
            continue

        # Skip empty or non-data lines
        if not line or line.startswith("FREQUENCY"):
            continue

        # Parse numeric rows
        try:
            parts = line.split()
            if len(parts) == 3:
                freq, mag, phase = map(float, parts)
                if current_section == "on":
                    chopper_on.append([freq, mag, phase])
                elif current_section == "off":
                    chopper_off.append([freq, mag, phase])
        except ValueError:
            continue

    # Convert to numpy arrays
    chopper_on = np.array(chopper_on)
    chopper_off = np.array(chopper_off)

    calibration_data['chopper_on'] = chopper_on
    calibration_data['chopper_off'] = chopper_off

    return calibration_data
