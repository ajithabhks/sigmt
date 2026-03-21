import re

import pandas as pd


def dms_to_decimal(dms_text: str) -> float:
    """
    Convert strings like:
      32 47'11.2"N
      077 16'53.5"E
    to decimal degrees.
    """
    m = re.match(r'\s*(\d+)\s+(\d+)\'([\d.]+)"([NSEW])\s*', dms_text)
    if not m:
        raise ValueError(f"Could not parse DMS value: {dms_text!r}")

    deg, minute, sec, hemi = m.groups()
    value = float(deg) + float(minute) / 60 + float(sec) / 3600

    if hemi in ("S", "W"):
        value = -value
    return value


def read_geomag_ascii(path: str):
    header = {}
    data_rows = []
    column_names = [
        "year", "month", "day",
        "hour", "minute", "second",
        "hx", "hy", "hz",
        "ex", "ey", "t_s", "t_e"
    ]

    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            if not stripped:
                continue

            # Header / metadata lines
            if stripped.startswith(";"):
                meta = stripped.lstrip(";").strip()
                if not meta:
                    continue

                # First line: "MS:GEOMAG-02 ..."
                if meta.startswith("MS:"):
                    header["MS"] = meta.split(":", 1)[1].strip()

                # Date + Time on same line
                elif meta.startswith("Date:"):
                    m = re.search(r"Date:\s*([^;]+)\s*;\s*Time:\s*(.+)$", meta)
                    if m:
                        header["file_date"] = m.group(1).strip()
                        header["file_time"] = m.group(2).strip()

                # Sampling
                elif meta.startswith("Sampling:"):
                    m = re.search(r"Sampling:\s*([\d.]+)", meta)
                    if m:
                        header["sampling_frequency"] = float(m.group(1))

                # Lat / Lon / Alt
                elif meta.startswith("Latitude:"):
                    m = re.search(
                        r"Latitude:\s*(.*?)\s*;\s*Longitude:\s*(.*?)\s*;\s*Altitude:\s*([+\-]?\d+(?:\.\d+)?)m",
                        meta
                    )
                    if m:
                        lat_text, lon_text, alt_text = m.groups()
                        header["latitude_deg"] = dms_to_decimal(lat_text)
                        header["longitude_deg"] = dms_to_decimal(lon_text)
                        header["altitude_m"] = float(alt_text)

                # Total field
                elif meta.startswith("Total Field:"):
                    m = re.search(
                        r"X\s*=\s*([+\-]?\d+)nT;\s*Y\s*=\s*([+\-]?\d+)nT;\s*Z\s*=\s*([+\-]?\d+)nT",
                        meta
                    )
                    if m:
                        header["total_field_x_nT"] = float(m.group(1))
                        header["total_field_y_nT"] = float(m.group(2))
                        header["total_field_z_nT"] = float(m.group(3))

                continue

            # Skip printed text header row like:
            # Date Time X [nT] ...
            if stripped.startswith("Date"):
                continue

            # Data rows begin with year
            if re.match(r"^\d{4}\s+\d{2}\s+\d{2}\s+\d{2}\s+\d{2}\s+[\d.]+", stripped):
                parts = re.split(r"\s+", stripped)
                if len(parts) >= 13:
                    data_rows.append(parts[:13])

    df = pd.DataFrame(data_rows, columns=column_names)

    # Convert numeric columns
    numeric_cols = column_names
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)

    # Build timestamp
    sec_int = df["second"].astype(int)
    sec_frac = df["second"] - sec_int

    df["timestamp"] = (
            pd.to_datetime(
                dict(
                    year=df["year"],
                    month=df["month"],
                    day=df["day"],
                    hour=df["hour"],
                    minute=df["minute"],
                    second=sec_int,
                )
            )
            + pd.to_timedelta(sec_frac, unit="s")
    )

    df["timestamp"] = df["timestamp"].dt.tz_localize("UTC")

    # Put timestamp first
    df = df[
        ["timestamp", "year", "month", "day", "hour", "minute", "second",
         "hx", "hy", "hz", "ex", "ey", "t_s", "t_e"]
    ]

    return header, df
