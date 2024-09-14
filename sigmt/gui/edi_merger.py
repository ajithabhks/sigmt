"""
GUI of EDI merger and does merging
"""

import numpy as np
import seaborn as sns
import xarray as xr
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QPushButton, QVBoxLayout, QWidget, QFileDialog
from matplotlib import pyplot as plt

from sigmt.utils.edi import edi_ops

if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)


class EDIMerger(QWidget):
    """
    GUI for EDI merger
    """

    def __init__(self):
        super().__init__()

        self.fig = None
        self.axs = None
        self.scatter1 = None
        self.scatter2 = None
        self.scatter3 = None
        self.scatter4 = None
        self.color1 = None
        self.initial_color1 = None
        self.initial_data = None
        self.scatter_list = None
        self.data = None
        self.edi_cont = None
        self.hmeas = None
        self.components = None
        self.edi_files = None
        self.setWindowTitle("EDI Merger")
        self.setGeometry(100, 100, 250, 200)
        self.setWindowIcon(QIcon(r'sigmt\images\sigmt.ico'))

        # Layout for buttons
        layout = QVBoxLayout()

        self.select_files_button = QPushButton("Select EDI Files")
        self.select_files_button.clicked.connect(self.select_files)
        layout.addWidget(self.select_files_button)

        self.plot_data_button = QPushButton("Plot Data")
        self.plot_data_button.clicked.connect(self.plot_data)
        layout.addWidget(self.plot_data_button)

        self.save_as_edi_button = QPushButton("Save as EDI")
        self.save_as_edi_button.clicked.connect(self.save_as_edi)
        layout.addWidget(self.save_as_edi_button)

        self.close_button = QPushButton("Close")
        self.close_button.clicked.connect(self.close_window)
        layout.addWidget(self.close_button)

        # Set the layout
        self.setLayout(layout)

    def select_files(self):
        """
        Open file dialog with .edi filter and allow multiple selections
        """
        self.components = [
            'freqs',
            'zxx_r', 'zxx_i', 'zxx_var',
            'zxy_r', 'zxy_i', 'zxy_var',
            'zyx_r', 'zyx_i', 'zyx_var',
            'zyy_r', 'zyy_i', 'zyy_var',
            'tzx_r', 'tzx_i', 'tzx_var',
            'tzy_r', 'tzy_i', 'tzy_var',
            'coh_ex', 'coh_ey'
        ]
        data = {}
        for component in self.components:
            data[component] = []
        #
        options = QFileDialog.Options()
        self.edi_files, _ = QFileDialog.getOpenFileNames(self, "Select EDI Files", "",
                                                         "EDI Files (*.edi);;"
                                                         "All Files (*)", options=options)
        if self.edi_files:
            colorx = self.generate_colors1()
            self.color1 = []
            self.hmeas = []
            cind = 0
            for edi_file in self.edi_files:
                self.edi_cont, data_temp, hmeas_temp = edi_ops.read_edi(edi_file)
                self.hmeas.append(hmeas_temp)
                for component in self.components:
                    data[component].extend(data_temp[component])
                for _ in range(len(data_temp['freqs'])):
                    self.color1.append(colorx[cind])
                cind = cind + 1
            # Rejects components with no data, but accepts if some are NaNs
            filtered_data = {key: value for key, value in data.items()
                             if key == 'freqs' or not np.all(np.isnan(value))}
            self.data = xr.Dataset(
                {key: (['frequency'], filtered_data[key]) for key in filtered_data},
                coords={'frequency': filtered_data['freqs']}
            )
            self.get_res_phase()
            if 5 in self.hmeas:
                indices_of_5 = [index for index, value in enumerate(self.hmeas) if value == 5]
                self.edi_cont, data_temp, hmeas_temp = edi_ops.read_edi(
                    self.edi_files[indices_of_5[0]])

    def plot_data(self):
        """
        Docs
        """
        self.fig, self.axs = plt.subplots(2)
        # Plot the scatter points
        self.scatter1 = self.axs[0].scatter(self.data['freqs'], self.data['rho_xy'],
                                            marker="s", color=self.color1,
                                            s=50,
                                            picker=True,
                                            label='XY')  # rho xy
        self.scatter2 = self.axs[0].scatter(self.data['freqs'], self.data['rho_yx'],
                                            marker="o", color=self.color1,
                                            s=50,
                                            picker=True,
                                            label='YX')  # rho yx
        self.scatter3 = self.axs[1].scatter(self.data['freqs'], self.data['phase_xy'],
                                            marker="s", color=self.color1,
                                            s=50,
                                            picker=True)  # phase xy
        self.scatter4 = self.axs[1].scatter(self.data['freqs'], self.data['phase_yx'],
                                            marker="o", color=self.color1,
                                            s=50,
                                            picker=True)  # phase yx

        # Store the scatter plots in a list for easy access
        self.scatter_list = [self.scatter1, self.scatter2, self.scatter3, self.scatter4]
        # Set log scale for scatter1 and scatter2
        self.axs[0].set_xscale('log')
        self.axs[0].set_yscale('log')

        # Add x-axis and y-axis labels
        self.axs[0].set_xlabel('Frequency (Hz)')
        self.axs[0].set_ylabel('App. Res. (Ohm.m.)')

        self.axs[1].set_xlabel('Frequency (Hz)')
        self.axs[1].set_ylabel('Phase (Deg.)')

        # Set log scale for x-axis of scatter3 and scatter4
        self.axs[1].set_xscale('log')
        self.axs[0].grid(which='both', linestyle='-.', linewidth=0.4)
        self.axs[1].grid(which='both', linestyle='-.', linewidth=0.4)
        self.axs[0].invert_xaxis()
        self.axs[1].invert_xaxis()
        self.axs[0].legend(loc='upper right', bbox_to_anchor=(1, 1))

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_right_click)

        plt.show(block=False)

    def on_pick(self, event):
        """
        Action for left click
        """
        # Get the index of the clicked point
        index = event.ind[0]

        # Determine which scatter plot was clicked
        for scatter in self.scatter_list:
            if event.artist == scatter:
                # Store initial state for undo
                self.initial_data = self.data.copy()
                self.initial_color1 = self.color1.copy()

                # Deleting...
                self.data = self.data.isel(
                    frequency=[i for i in range(len(self.data['frequency'])) if i != index])
                del self.color1[index]

                # Update the scatter plots with the new data
                self.scatter1.set_offsets(
                    np.column_stack((self.data['freqs'], self.data['rho_xy'])))
                self.scatter2.set_offsets(
                    np.column_stack((self.data['freqs'], self.data['rho_yx'])))
                self.scatter3.set_offsets(
                    np.column_stack((self.data['freqs'], self.data['phase_xy'])))
                self.scatter4.set_offsets(
                    np.column_stack((self.data['freqs'], self.data['phase_yx'])))

                # Update the colors
                self.scatter1.set_color(self.color1)
                self.scatter2.set_color(self.color1)
                self.scatter3.set_color(self.color1)
                self.scatter4.set_color(self.color1)

                # Update the plots
                self.axs[0].figure.canvas.draw()

    def on_right_click(self, event):
        """
        Action for right click
        """
        if event.button == 3:  # Right mouse button click
            self.data = self.initial_data.copy()
            self.color1 = self.initial_color1.copy()
            # Restore original scatter plots
            self.scatter1.set_offsets(np.column_stack((self.data['freqs'], self.data['rho_xy'])))
            self.scatter2.set_offsets(np.column_stack((self.data['freqs'], self.data['rho_yx'])))
            self.scatter3.set_offsets(np.column_stack((self.data['freqs'], self.data['phase_xy'])))
            self.scatter4.set_offsets(np.column_stack((self.data['freqs'], self.data['phase_yx'])))

            self.scatter1.set_color(self.color1)
            self.scatter2.set_color(self.color1)
            self.scatter3.set_color(self.color1)
            self.scatter4.set_color(self.color1)
            # Update the plots
            self.axs[0].figure.canvas.draw()

    def save_as_edi(self):
        """
        Docs
        """
        data_frame = self.data.to_dataframe().reset_index()
        data_frame = data_frame.sort_values(by=["frequency"], ascending=[False])

        file_formats = "EDI Files (*.edi);;All Files (*)"
        save_path, _ = QFileDialog.getSaveFileName(None, "Save File", "", file_formats)
        if save_path:  # If user selects a path
            with open(save_path, 'w', encoding='utf-8') as file_handle:
                nfindex = 0
                for j in self.edi_cont:
                    if j[:8] == '  NFREQ=':
                        break
                    nfindex += 1
                for i in range(nfindex):
                    file_handle.write(self.edi_cont[i])
                file_handle.write("  NFREQ=" + str(np.shape(self.data['freqs'])[0]) + "\n")
                findex = 0
                for j in self.edi_cont:
                    if j[:5] == '>FREQ':
                        break
                    findex += 1
                for i in range(nfindex + 1, findex - 1):
                    file_handle.write(self.edi_cont[i])

                # ===== FREQ section =====
                file_handle.write("\n")
                file_handle.write(">FREQ //" + str(len(data_frame['frequency'])))
                file_handle.write("\n")

                for i in range(len(data_frame['frequency'])):
                    k = i + 1
                    file_handle.write(f" {data_frame['frequency'][i]:.9E} ")
                    if (k % 6 == 0) and (k != 0):
                        file_handle.write("\n")

                if any(col in data_frame.columns for col in ['zxy_r', 'zyx_r']):
                    # ===== ZROT section =====
                    file_handle.write("\n \n")
                    file_handle.write(">ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        file_handle.write(f" {0:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXXR section =====
                if 'zxx_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXXR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxx_r'][i] < 0:
                            file_handle.write(f"{data_frame['zxx_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxx_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXXI section =====
                if 'zxx_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXXI ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxx_i'][i] < 0:
                            file_handle.write(f"{data_frame['zxx_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxx_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXX.VAR section =====
                if 'zxx_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXX.VAR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxx_var'][i] < 0:
                            file_handle.write(f"{data_frame['zxx_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxx_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXYR section =====
                if 'zxy_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXYR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxy_r'][i] < 0:
                            file_handle.write(f"{data_frame['zxy_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxy_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXYI section =====
                if 'zxy_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXYI ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxy_i'][i] < 0:
                            file_handle.write(f"{data_frame['zxy_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxy_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZXY.VAR section =====
                if 'zxy_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZXY.VAR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zxy_var'][i] < 0:
                            file_handle.write(f"{data_frame['zxy_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zxy_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYXR section =====
                if 'zyx_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYXR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyx_r'][i] < 0:
                            file_handle.write(f"{data_frame['zyx_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyx_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYXI section =====
                if 'zyx_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYXI ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyx_i'][i] < 0:
                            file_handle.write(f"{data_frame['zyx_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyx_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYX.VAR section =====
                if 'zyx_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYX.VAR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyx_var'][i] < 0:
                            file_handle.write(f"{data_frame['zyx_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyx_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYYR section =====
                if 'zyy_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYYR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyy_r'][i] < 0:
                            file_handle.write(f"{data_frame['zyy_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyy_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYYI section =====
                if 'zyy_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYYI ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyy_i'][i] < 0:
                            file_handle.write(f"{data_frame['zyy_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyy_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== ZYY.VAR section =====
                if 'zyy_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">ZYY.VAR ROT=ZROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['zyy_var'][i] < 0:
                            file_handle.write(f"{data_frame['zyy_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['zyy_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TROT section =====
                if any(col in data_frame.columns for col in ['tzx_r', 'tzy_r']):
                    file_handle.write("\n \n")
                    file_handle.write(">TROT.EXP //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        file_handle.write(f" {0:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TXR.EXP section =====
                if 'tzx_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TXR.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzx_r'][i] < 0:
                            file_handle.write(f"{data_frame['tzx_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzx_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TXI.EXP section =====
                if 'tzx_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TXI.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzx_i'][i] < 0:
                            file_handle.write(f"{data_frame['tzx_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzx_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TXVAR.EXP section =====
                if 'tzx_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TXVAR.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzx_var'][i] < 0:
                            file_handle.write(f"{data_frame['tzx_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzx_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TYR.EXP section =====
                if 'tzy_r' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TYR.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzy_r'][i] < 0:
                            file_handle.write(f"{data_frame['tzy_r'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzy_r'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TYI.EXP section =====
                if 'tzy_i' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TYI.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzy_i'][i] < 0:
                            file_handle.write(f"{data_frame['tzy_i'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzy_i'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== TYVAR.EXP section =====
                if 'tzy_var' in data_frame.columns:
                    file_handle.write("\n \n")
                    file_handle.write(">TYVAR.EXP ROT=TROT //" + str(len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        if data_frame['tzy_var'][i] < 0:
                            file_handle.write(f"{data_frame['tzy_var'][i]:.9E} ")
                        else:
                            file_handle.write(f" {data_frame['tzy_var'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== cohEx section =====
                if 'coh_ex' in data_frame.columns:
                    file_handle.write("\n \n")
                    if 5 in self.hmeas:
                        file_handle.write(
                            ">COH  MEAS1=" + str(101.001) +
                            " MEAS2=" + str(107.001) + "  ROT=NORTH //" + str(
                                len(data_frame['frequency'])))
                    else:
                        file_handle.write(
                            ">COH  MEAS1=" + str(101.001) +
                            " MEAS2=" + str(105.001) + "  ROT=NORTH //" + str(
                                len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        file_handle.write(f" {data_frame['coh_ex'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")

                # ===== cohEy section =====
                if 'coh_ey' in data_frame.columns:
                    file_handle.write("\n \n")
                    if 5 in self.hmeas:
                        file_handle.write(
                            ">COH  MEAS1=" + str(102.001) +
                            " MEAS2=" + str(106.001) + "  ROT=NORTH //" + str(
                                len(data_frame['frequency'])))
                    else:
                        file_handle.write(
                            ">COH  MEAS1=" + str(102.001) +
                            " MEAS2=" + str(104.001) + "  ROT=NORTH //" + str(
                                len(data_frame['frequency'])))
                    file_handle.write("\n")

                    for i in range(len(data_frame['frequency'])):
                        k = i + 1
                        file_handle.write(f" {data_frame['coh_ey'][i]:.9E} ")
                        if (k % 6 == 0) and (k != 0):
                            file_handle.write("\n")
                    file_handle.write("\n \n")
                    file_handle.write(">END")
                    file_handle.close()

    def generate_colors1(self):
        """
        Generate a list of 'num_colors' distinct colors
        """
        return sns.color_palette("dark", len(self.edi_files))

    def generate_colors2(self):
        """
        Generate a list of 'num_colors' distinct colors
        """
        return sns.color_palette("deep", len(self.edi_files))

    def get_res_phase(self):
        """
        Calculates resistivity and phase from impedance values
        """
        ftlist = self.data['freqs']
        zxy = self.data['zxy_r'] + 1j * self.data['zxy_i']
        zyx = self.data['zyx_r'] + 1j * self.data['zyx_i']

        self.data['rho_xy'] = (0.2 / ftlist) * (abs(zxy) ** 2)
        self.data['rho_yx'] = (0.2 / ftlist) * (abs(zyx) ** 2)
        self.data['phase_xy'] = np.degrees(np.arctan(zxy.imag / zxy.real))
        self.data['phase_yx'] = np.degrees(np.arctan(zyx.imag / zyx.real))

    def close_window(self):
        """
        Closes the window
        """
        self.close()
