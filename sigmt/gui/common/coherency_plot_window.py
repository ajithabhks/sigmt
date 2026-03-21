import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


def _get_custom_colormap():
    return plt.cm.viridis


def _get_tf_key(coh_keys, index):
    mapping = {
        "coh_ex": "zxx",
        "coh_ey": "zxy",
        "coh_hz": "tzx",
    }
    return mapping[coh_keys[index]]


class CoherencyPlotWindow(QWidget):
    def __init__(self, dataset: xr.Dataset, parent=None):
        super().__init__(parent)

        if dataset is None or "frequency" not in dataset:
            raise ValueError("dataset invalid or missing 'frequency'")

        self.dataset = dataset
        self.custom_colormap = _get_custom_colormap()

        self.coh_keys = []
        if "coh_ex" in dataset:
            self.coh_keys.append("coh_ex")
        if "coh_ey" in dataset:
            self.coh_keys.append("coh_ey")
        if "coh_hz" in dataset:
            self.coh_keys.append("coh_hz")

        if not self.coh_keys:
            raise ValueError("No coherency variables found in dataset")

        self.setWindowTitle("Coherency Viewer")
        self.resize(1400, 800)

        self.freq_label = QLabel("Frequency:")
        self.freq_combo = QComboBox()

        self.frequencies = self.dataset["frequency"].values
        for freq in self.frequencies:
            self.freq_combo.addItem(f"{freq:.6f} Hz", freq)

        self.freq_combo.currentIndexChanged.connect(self.update_plot)

        top_layout = QHBoxLayout()
        top_layout.addWidget(self.freq_label)
        top_layout.addWidget(self.freq_combo)
        top_layout.addStretch()

        self.figure = Figure(constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.canvas)
        self.setLayout(main_layout)

        self.update_plot()

    @staticmethod
    def _clean_xy(x, y):
        x = np.asarray(x).ravel()
        y = np.asarray(y).ravel()
        mask = np.isfinite(x) & np.isfinite(y)
        return x[mask], y[mask]

    @staticmethod
    def _set_square_equal(ax, x, y, pad_frac=0.05):
        x, y = CoherencyPlotWindow._clean_xy(x, y)

        if x.size == 0:
            ax.set_box_aspect(1)
            ax.set_aspect("equal", adjustable="datalim")
            return

        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()

        xmid = 0.5 * (xmin + xmax)
        ymid = 0.5 * (ymin + ymax)

        half = max(xmax - xmin, ymax - ymin) / 2.0
        if half == 0:
            half = 1.0

        half *= (1.0 + pad_frac)

        ax.set_xlim(xmid - half, xmid + half)
        ax.set_ylim(ymid - half, ymid + half)

        # Square drawing area + equal data scaling
        ax.set_box_aspect(1)
        ax.set_aspect("equal", adjustable="box")

    def update_plot(self):
        freq = self.freq_combo.currentData()
        data_for_freq = self.dataset.sel(frequency=freq)

        self.figure.clear()
        num_coh = len(self.coh_keys)

        # generous spacing handled by constrained_layout
        axes = self.figure.subplots(1, num_coh, squeeze=False)[0]

        for i, ax in enumerate(axes):
            tf_component = _get_tf_key(self.coh_keys, i)

            x = data_for_freq[tf_component].real.values
            y = data_for_freq[tf_component].imag.values
            c = data_for_freq[self.coh_keys[i]].values

            sc = ax.scatter(
                x,
                y,
                c=c,
                cmap=self.custom_colormap,
                vmin=0,
                vmax=1,
                s=35,
                edgecolors="none",
            )

            ax.set_title(self.coh_keys[i], fontsize=14)
            ax.set_xlabel(f"Real ({tf_component})", fontsize=11)
            ax.set_ylabel(f"Imag ({tf_component})", fontsize=11)

            self._set_square_equal(ax, x, y)

            # slimmer colorbar so it doesn't eat the plot area
            cbar = self.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(self.coh_keys[i], fontsize=11)

            ax.grid(True, alpha=0.25)

        self.figure.suptitle(f"Coherency plots for {freq:.6f} Hz", fontsize=16)
        self.canvas.draw()
