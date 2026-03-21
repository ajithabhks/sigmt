import numpy as np
import xarray as xr

from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure


class PolarizationDirectionWindow(QWidget):
    def __init__(self, dataset: xr.Dataset, parent=None):
        super().__init__(parent)

        if dataset is None or "frequency" not in dataset:
            raise ValueError("dataset invalid or missing 'frequency'")

        if "alpha_h" not in dataset or "alpha_e" not in dataset:
            raise ValueError("dataset must contain 'alpha_h' and 'alpha_e'")

        self.dataset = dataset

        self.setWindowTitle("Polarization Direction Viewer")
        self.resize(1100, 800)

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
        self.toolbar = NavigationToolbar(self.canvas, self)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)
        self.setLayout(main_layout)

        self.update_plot()

    def update_plot(self):
        freq = self.freq_combo.currentData()
        data_for_freq = self.dataset.sel(frequency=freq)

        self.figure.clear()

        ax1, ax2 = self.figure.subplots(2, 1, sharex=True)

        time_array = np.arange(data_for_freq.coords["time_window"].size)
        yticks = np.arange(-90, 100, 20)

        alpha_h = np.asarray(data_for_freq["alpha_h"].values).ravel()
        alpha_e = np.asarray(data_for_freq["alpha_e"].values).ravel()

        ax1.scatter(time_array, alpha_h, s=20)
        ax1.set_ylim(-90, 90)
        ax1.set_ylabel(r"$\alpha_H$ (deg.)")
        ax1.set_yticks(yticks)
        ax1.grid(True, alpha=0.25)

        ax2.scatter(time_array, alpha_e, s=20)
        ax2.set_ylim(-90, 90)
        ax2.set_ylabel(r"$\alpha_E$ (deg.)")
        ax2.set_xlabel("Time window")
        ax2.set_yticks(yticks)
        ax2.grid(True, alpha=0.25)

        self.canvas.draw()
