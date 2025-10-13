import struct
import numpy as np
from pathlib import Path


class PhoenixMTReader:
    """
    Reader for Phoenix Geophysics MTU/RXU/MTU-5D format files (v2.0+)
    """

    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.header = {}
        self.data = None

    def read_header(self):
        """
        Reads and parses the 128-byte header
        """
        with open(self.file_path, 'rb') as f:
            header_bytes = f.read(128)

        if len(header_bytes) < 128:
            raise ValueError("Header too short (expected 128 bytes)")

        file_type = header_bytes[0]
        file_version = header_bytes[1]
        header_len = struct.unpack('<H', header_bytes[2:4])[0]

        # Common fields for both continuous (type 1) and decimated (type 2)
        instrument = header_bytes[4:12].decode('ascii', errors='ignore').strip('\x00 ')
        serial = header_bytes[12:20].decode('ascii', errors='ignore').strip('\x00 ')
        record_id = struct.unpack('<I', header_bytes[20:24])[0]
        channel_id = header_bytes[24]
        file_seq = struct.unpack('<I', header_bytes[25:29])[0]
        sample_rate_base = struct.unpack('<H', header_bytes[59:61])[0]
        bytes_per_sample = header_bytes[62]

        self.header = {
            "file_type": file_type,
            "file_version": file_version,
            "header_length": header_len,
            "instrument": instrument,
            "serial_number": serial,
            "recording_id": record_id,
            "channel_id": channel_id,
            "file_sequence": file_seq,
            "sample_rate": sample_rate_base,
            "bytes_per_sample": bytes_per_sample,
        }

    def read_data(self):
        """
        Reads data payload after header
        """
        with open(self.file_path, 'rb') as f:
            f.seek(128)
            raw = f.read()

        if self.header["file_type"] == 1:
            # Continuous data: 64-byte frames = 20 * 3-byte samples + 4-byte footer
            frames = len(raw) // 64
            samples = np.empty(frames * 20, dtype=np.int32)
            for i in range(frames):
                frame = raw[i * 64:(i + 1) * 64]
                for j in range(20):
                    b = frame[j * 3:(j + 1) * 3]
                    val = int.from_bytes(b, byteorder='big', signed=True)
                    samples[i * 20 + j] = val
            self.data = samples

        elif self.header["file_type"] == 2:
            # Decimated data: float32 samples in Volts
            self.data = np.frombuffer(raw, dtype='<f4')

        else:
            raise ValueError(f"Unknown file_type {self.header['file_type']}")

    def read(self):
        """
        Wrapper
        """
        self.read_header()
        self.read_data()
        return self.header, self.data


# reader = PhoenixMTReader("AAAAA_BBBBBBBB_C_DDDDDDDD.td_24k")
# header, data = reader.read()
