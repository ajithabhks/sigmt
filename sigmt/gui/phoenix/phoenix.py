"""
First landing window for the Phoenix specific operations.

"""

import json
import os
import re
import time

import numpy as np
import xarray as xr
import yaml
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (
    QMainWindow, QAction, QFileDialog, QMessageBox,
    QWidget, QVBoxLayout, QComboBox, QLabel,
    QPushButton, QHBoxLayout, QRadioButton, QGroupBox,
    QGridLayout, QLineEdit, QSizePolicy, QApplication, QProgressDialog
)
from scipy import signal

from sigmt.core import data_selection_tools as dstools
from sigmt.core import perform_data_selection as pds
from sigmt.core import plots
from sigmt.core.band_averaging import BandAveraging
from sigmt.core.robust_estimation import RobustEstimation
from sigmt.gui.about_dialog import AboutDialog
from sigmt.gui.disclaimer_dialog import DisclaimerDialog
from sigmt.gui.edi_merger import EDIMerger
from sigmt.gui.project_related.create_project import ProjectSetupDialog
from sigmt.gui.project_related.edit_project import EditProjectSetupDialog
from sigmt.utils import utils
from sigmt.utils.edi import edi_ops
from sigmt.utils.phoenix import phoenix_utils

if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)


class MainWindow(QMainWindow):
    """
    Definition of the main window for Metronix processing

    """

    def __init__(self):
        """
        Constructor

        """
        super().__init__()
        self.notch_status = None
        self.md_threshold_entry = None
        self.coh_plot_button = None
        self.pd_plot_button = None
        self.about_dialog = None
        self.read_time_series_button = None
        self.decimation_list_dropdown = None
        self.decimate_button = None
        self.new_fs = None
        self.notch_radio_off = None
        self.notch_radio_on = None
        self.fft_length_dropdown = None
        self.parzen_radius_entry = None
        self.save_parameters_button = None
        self.perform_bandavg_button = None
        self.coherency_threshold_edit = None
        self.minimum_data_edit = None
        self.apply_coh_thresh_button = None
        self.clear_coh_thresh_button = None
        self.pd_min_edit = None
        self.pd_max_edit = None
        self.pd_combo_box = None
        self.apply_pd_thresh_button = None
        self.clear_pd_thresh_button = None
        self.plot_impedance_button = None
        self.plot_tipper_button = None
        self.plot_coherency_button = None
        self.save_as_edi_button = None
        self.open_edi_merger_button = None
        self.allsites = None
        self.localsite_path = None
        self.remotesite_path = None
        self.sfreq_selected = None
        self.menubar = None
        self.edi_merger = None
        self.plot_edi_button = None
        self.estimates = None
        self.perform_robust_estimation_button = None
        self.bandavg_dataset = None
        self.procinfo = {}
        self.project_dir = None
        self.project_setup = None
        self.localsite = None
        self.localsite_meas = None
        self.localsite_dropdown = None
        self.remotesite_dropdown = None
        self.sampling_frequency_dropdown = None
        self.disclaimer_dialog = None
        self.fft_values = ['262144', '131072', '65536', '32768', '16384',
                           '8192', '4096', '2048', '1024', '512', '256']
        #

        self.time_series = None
        self.calibration_data_electric = None
        self.calibration_data_magnetic = None

        self.recmeta_data_local = None
        self.recmeta_data_remote = None
        self.file_extension = None
        self.file_type = None

        self.remotesite = None
        #
        self.interface = 'Phoenix'
        self.init_ui()

    def init_ui(self) -> None:
        """
        Initial display user interface

        :return: None
        :rtype: NoneType

        """
        self.setWindowTitle('[No Project Opened] SigMT | A Tool for '
                            f'Magnetotelluric Data Processing ({self.interface})')
        self.setWindowIcon(QIcon(r'sigmt\images\sigmt.ico'))
        self.setGeometry(100, 100, 700, 500)

        # Creating a menubar
        self.menubar = self.menuBar()

        # Creating File menubar item
        # -------------------------------
        file = self.menubar.addMenu("File")

        create_project_action = QAction("Create a Project", self)
        create_project_action.triggered.connect(self.create_project)
        file.addAction(create_project_action)

        open_project_action = QAction("Open a Project", self)
        open_project_action.triggered.connect(self.open_project)
        file.addAction(open_project_action)

        close_project_action = QAction("Close Project", self)
        close_project_action.triggered.connect(self.close_project)
        file.addAction(close_project_action)

        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.exit_all)
        file.addAction(exit_action)

        # Creating Edit menubar item
        # -------------------------------
        edit = self.menubar.addMenu("Edit")

        edit_project_action = QAction("Edit Project Setup", self)
        edit_project_action.triggered.connect(self.edit_project_setup)
        edit.addAction(edit_project_action)

        about_menu = self.menubar.addMenu("About")
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about_dialog)
        about_menu.addAction(about_action)

        disclaimer_menu = self.menubar.addMenu("⚠️ IMPORTANT DISCLAIMER")
        disclaimer_menu.setProperty("warning", True)

        disclaimer_menu.setStyleSheet("""
            QMenu {
                background-color: #fff5f5;
                color: #842029;
                font-weight: bold;
            }
            QMenu::item:selected {
                background-color: #f1aeb5;
            }
        """)

        disclaimer_action = QAction("Read this BEFORE using Phoenix", self)
        disclaimer_action.triggered.connect(self.show_disclaimer)
        disclaimer_menu.addAction(disclaimer_action)

        # Central Widget
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QVBoxLayout()

        # Subsection 1
        section1 = QGroupBox("Selection of sites and sampling frequency")
        section1_layout = QVBoxLayout()
        section1.setLayout(section1_layout)
        #
        load_sites_widget = QWidget()
        load_sites_layout = QHBoxLayout()
        load_sites_layout.setContentsMargins(0, 0, 0, 0)
        load_sites_widget.setLayout(load_sites_layout)
        load_sites_button = QPushButton("Load Sites")
        load_sites_button.clicked.connect(self.load_sites)
        load_sites_layout.addWidget(load_sites_button)
        load_sites_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        section1_layout.addWidget(load_sites_widget)
        #
        site_selection_widget = QWidget()
        site_selection_layout = QHBoxLayout()
        site_selection_layout.setContentsMargins(0, 0, 0, 0)
        site_selection_widget.setLayout(site_selection_layout)
        site_selection_layout.addWidget(QLabel("Select Local Site:"))
        self.localsite_dropdown = QComboBox()
        self.localsite_dropdown.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        site_selection_layout.addWidget(self.localsite_dropdown)
        site_selection_layout.addWidget(QLabel("Select Remote Site:"))
        self.remotesite_dropdown = QComboBox()
        self.remotesite_dropdown.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        site_selection_layout.addWidget(self.remotesite_dropdown)
        section1_layout.addWidget(site_selection_widget)
        #
        section1_grid_widget = QWidget()
        section1_grid_layout = QGridLayout()
        section1_grid_layout.setContentsMargins(0, 0, 0, 0)
        section1_grid_widget.setLayout(section1_grid_layout)
        section1_grid_layout.addWidget(QLabel("Select Sampling Frequency (Hz)"), 0, 0)
        self.sampling_frequency_dropdown = QComboBox()
        section1_grid_layout.addWidget(self.sampling_frequency_dropdown, 0, 1, 1, 40)
        self.read_time_series_button = QPushButton("Read time series")
        self.read_time_series_button.clicked.connect(self.read_ts)
        section1_grid_layout.addWidget(self.read_time_series_button, 1, 0)
        section1_grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        #
        section1_layout.addWidget(section1_grid_widget)
        section1_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Subsection 2
        section2 = QGroupBox("Decimation")
        section2_layout = QGridLayout()
        section2_layout.addWidget(QLabel("Select Decimation from the list:"), 0, 0)
        self.decimation_list_dropdown = QComboBox()
        self.decimation_list_dropdown.addItems(['5'])
        section2_layout.addWidget(self.decimation_list_dropdown, 0, 1)
        self.decimate_button = QPushButton("Decimate")
        self.decimate_button.clicked.connect(self.decimate)
        section2_layout.addWidget(self.decimate_button, 0, 2)
        self.new_fs = QLabel()
        section2_layout.addWidget(self.new_fs, 0, 3)
        section2.setLayout(section2_layout)
        section2_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Subsection 3
        section3 = QGroupBox("Parameter Selection")
        section3_layout = QGridLayout()
        section3_layout.addWidget(QLabel("Notch filter:"), 0, 0)
        # ----Notch radio button layout starts----
        self.notch_status = 'off'
        notch_radio_layout = QHBoxLayout()
        self.notch_radio_off = QRadioButton("Off")
        self.notch_radio_on = QRadioButton("On")
        self.notch_radio_off.setChecked(True)
        notch_radio_layout.addWidget(self.notch_radio_off)
        notch_radio_layout.addWidget(self.notch_radio_on)
        #
        self.notch_radio_off.clicked.connect(self.handle_radio_button)
        self.notch_radio_on.clicked.connect(self.handle_radio_button)
        # ----Notch radio button layout ends----
        section3_layout.addLayout(notch_radio_layout, 0, 1)
        section3_layout.addWidget(QLabel("FFT Length"), 1, 0)
        self.fft_length_dropdown = QComboBox()
        section3_layout.addWidget(self.fft_length_dropdown, 1, 1)
        section3_layout.addWidget(QLabel("Parzen window radius"), 1, 2)
        self.parzen_radius_entry = QLineEdit()
        section3_layout.addWidget(self.parzen_radius_entry, 1, 3)
        section3_layout.addWidget(QLabel("Mahalanobis Distance threshold"), 2, 0)
        self.md_threshold_entry = QLineEdit()
        section3_layout.addWidget(self.md_threshold_entry, 2, 1)
        self.save_parameters_button = QPushButton("Save Parameters")
        self.save_parameters_button.clicked.connect(self.save_parameters)
        section3_layout.addWidget(self.save_parameters_button, 3, 0)
        #
        section3.setLayout(section3_layout)
        section3_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Subsection 4
        section4 = QGroupBox("Band Averaging")
        section4_layout = QGridLayout()
        self.perform_bandavg_button = QPushButton("Perform Band Averaging")
        self.perform_bandavg_button.clicked.connect(self.perform_bandavg)
        section4_layout.addWidget(self.perform_bandavg_button, 0, 0)
        section4.setLayout(section4_layout)
        section4_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Subsection 5
        section5 = QGroupBox("Data Selection")
        section5_layout = QVBoxLayout()
        # --- Coherency ---
        section5_coh_layout = QHBoxLayout()
        self.coh_plot_button = QPushButton("Check Coherency Plots")
        self.coh_plot_button.clicked.connect(self.plot_coh_all)
        section5_coh_layout.addWidget(self.coh_plot_button)
        section5_coh_layout.addWidget(QLabel("Coherency Threshold"))
        self.coherency_threshold_edit = QLineEdit()
        self.coherency_threshold_edit.setText("0.9")
        section5_coh_layout.addWidget(self.coherency_threshold_edit)
        section5_coh_layout.addWidget(QLabel("Min % Data"))
        self.minimum_data_edit = QLineEdit()
        self.minimum_data_edit.setText("20")
        section5_coh_layout.addWidget(self.minimum_data_edit)
        self.apply_coh_thresh_button = QPushButton("Apply coherency threshold")
        self.apply_coh_thresh_button.clicked.connect(self.apply_coh_thresh)
        section5_coh_layout.addWidget(self.apply_coh_thresh_button)
        self.clear_coh_thresh_button = QPushButton("Clear")
        self.clear_coh_thresh_button.clicked.connect(self.clear_coh_thresh)
        section5_coh_layout.addWidget(self.clear_coh_thresh_button)
        section5_layout.addLayout(section5_coh_layout)
        # --- Polarization direction ---
        section5_pd_layout = QHBoxLayout()
        self.pd_plot_button = QPushButton("Check Polarization Direction (PD) Plots")
        self.pd_plot_button.clicked.connect(self.plot_pd_all)
        section5_pd_layout.addWidget(self.pd_plot_button)
        section5_pd_layout.addWidget(QLabel("PD Min:"))
        self.pd_min_edit = QLineEdit()
        self.pd_min_edit.setText("-10")
        section5_pd_layout.addWidget(self.pd_min_edit)
        section5_pd_layout.addWidget(QLabel("PD Max:"))
        self.pd_max_edit = QLineEdit()
        self.pd_max_edit.setText("10")
        section5_pd_layout.addWidget(self.pd_max_edit)
        section5_pd_layout.addWidget(QLabel("Component:"))
        self.pd_combo_box = QComboBox()
        self.pd_combo_box.addItem("Electric")
        self.pd_combo_box.addItem("Magnetic")
        section5_pd_layout.addWidget(self.pd_combo_box)
        self.apply_pd_thresh_button = QPushButton("Perform PD thresholding")
        self.apply_pd_thresh_button.clicked.connect(self.apply_pd_thresh)
        section5_pd_layout.addWidget(self.apply_pd_thresh_button)
        self.clear_pd_thresh_button = QPushButton("Clear")
        self.clear_pd_thresh_button.clicked.connect(self.clear_pd_thresh)
        section5_pd_layout.addWidget(self.clear_pd_thresh_button)
        section5_layout.addLayout(section5_pd_layout)
        # --------
        section5.setLayout(section5_layout)
        section5_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        club_six_seven_eight = QHBoxLayout()

        # Subsection 6
        section6 = QGroupBox("Robust Estimation")
        section6_layout = QHBoxLayout()
        self.perform_robust_estimation_button = QPushButton("Perform Robust Estimation")
        self.perform_robust_estimation_button.clicked.connect(self.perform_robust_estimation)
        section6_layout.addWidget(self.perform_robust_estimation_button)
        section6.setLayout(section6_layout)
        section6_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        club_six_seven_eight.addWidget(section6)

        # Subsection 7
        section7 = QGroupBox("Plots")
        section7_layout = QHBoxLayout()
        self.plot_impedance_button = QPushButton("App. Res and Phase")
        self.plot_impedance_button.clicked.connect(self.plot_impedance)
        section7_layout.addWidget(self.plot_impedance_button)
        self.plot_tipper_button = QPushButton("Tipper")
        section7_layout.addWidget(self.plot_tipper_button)
        self.plot_tipper_button.clicked.connect(self.plot_tipper)

        self.plot_coherency_button = QPushButton("Coherency")
        section7_layout.addWidget(self.plot_coherency_button)
        self.plot_coherency_button.clicked.connect(self.plot_coherency)

        section7.setLayout(section7_layout)
        section7_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        club_six_seven_eight.addWidget(section7)

        # Subsection 8
        section8 = QGroupBox("EDI Operations")
        section8_layout = QHBoxLayout()
        self.save_as_edi_button = QPushButton("Save as EDI")
        section8_layout.addWidget(self.save_as_edi_button)
        self.save_as_edi_button.clicked.connect(self.save_as_edi)
        #
        self.open_edi_merger_button = QPushButton("Open EDI merger")
        section8_layout.addWidget(self.open_edi_merger_button)
        self.open_edi_merger_button.clicked.connect(self.open_edi_merger)
        #
        self.plot_edi_button = QPushButton("Plot an EDI")
        section8_layout.addWidget(self.plot_edi_button)
        self.plot_edi_button.setEnabled(False)
        section8.setLayout(section8_layout)
        section8_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        club_six_seven_eight.addWidget(section8)

        club_six_seven_eight.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Adding subsections to the main layout
        main_layout.addWidget(section1)
        main_layout.addWidget(section2)
        main_layout.addWidget(section3)
        main_layout.addWidget(section4)
        main_layout.addWidget(section5)
        main_layout.addLayout(club_six_seven_eight)
        main_layout.setAlignment(Qt.AlignTop)

        # Set main layout to central widget
        central_widget.setLayout(main_layout)

    def create_project(self) -> None:
        """
        It opens the dialog box when the 'create project' option is
        selected from the menu bar.
        A project is mandatory for further processing.

        :return: None
        :rtype: NoneType

        """
        # Prompt user to select project directory
        self.project_dir = QFileDialog.getExistingDirectory(self, "Select Project Directory", "/")
        if self.project_dir:
            # Ask questions and save project setup
            dialog = ProjectSetupDialog(self, interface=self.interface)

            if dialog.exec_():
                # Create the necessary folders
                folders = ['time_series', 'edi', 'project_setup', 'calibration_files']
                for folder in folders:
                    folder_path = os.path.join(self.project_dir, folder)
                    os.makedirs(folder_path, exist_ok=True)

                self.project_setup = dialog.get_values()

                # Save setup data to YAML file
                setup_file_path = os.path.join(self.project_dir, 'project_setup', 'setup.yaml')
                with open(setup_file_path, 'w') as yaml_file:
                    yaml.dump(self.project_setup, yaml_file)
                QMessageBox.information(self, 'Done',
                                        f'Project created at: {self.project_dir}')
                QMessageBox.information(self, 'Information',
                                        f'Opening new project: {self.project_dir}')
                self.setWindowTitle(self.project_setup[
                                        'project_name'] +
                                    ' - SigMT | '
                                    f'A Tool for Magnetotelluric Data Processing ({self.interface})')

    def open_project(self) -> None:
        """
        To open an existing project.
        
        :return: None
        :rtype: NoneType

        """
        self.new_fs.hide()
        self.fft_length_dropdown.clear()
        self.parzen_radius_entry.clear()
        self.md_threshold_entry.clear()
        #
        # Prompt user to select project directory
        self.project_dir = QFileDialog.getExistingDirectory(
            self, "Select Project Directory", "/")
        try:
            self.project_setup = utils.read_yaml_file(os.path.join(
                self.project_dir, 'project_setup', 'setup.yaml'))
            # To support older projects.
            if 'target_frequency_table_type' not in self.project_setup:
                self.project_setup['target_frequency_table_type'] = 'Default'
            if self.project_setup['interface'] == self.interface:
                self.setWindowTitle(
                    self.project_setup[
                        'project_name'] +
                    ' - SigMT | A Tool for Magnetotelluric Data Processing '
                    f'({self.interface})'
                )
            else:
                QMessageBox.warning(self, "Warning", f"Not a {self.interface} Project")
        except:
            QMessageBox.warning(self, "Warning", "Not a valid SigMT Project")

    def close_project(self) -> None:
        """
        Closes the current project and re-opens the MainWindow.

        :return: None
        :rtype: NoneType

        """
        self.close()
        self.__init__()
        self.show()

    def edit_project_setup(self) -> None:
        """
        If user wants to edit the project setup, it can be done using the option in menu bar. This 
        method manages the project setup editing.

        :return: None
        :rtype: NoneType

        """
        if self.project_setup is not None:
            try:
                dialog = EditProjectSetupDialog(
                    self, interface=self.interface, project_setup=self.project_setup)
                if dialog.exec_():
                    self.project_setup = dialog.get_values()

                    # Save setup data to YAML file
                    setup_file_path = os.path.join(self.project_dir, 'project_setup', 'setup.yaml')
                    with open(setup_file_path, 'w') as yaml_file:
                        yaml.dump(self.project_setup, yaml_file)
                    QMessageBox.information(
                        self, 'Done', f'Project setup edited at: {self.project_dir}')
                self.setWindowTitle(
                    self.project_setup[
                        'project_name'] +
                    ' - SigMT | A Tool for Magnetotelluric Data Processing '
                    f'({self.interface})')
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred while editing setup: {e}")
        else:
            QMessageBox.warning(self, "Warning", "Please create and/or open a project!")

    def show_about_dialog(self):
        """
        Opens a dialog box with some information about SigMT and version.

        :return: None
        :rtype: NoneType

        """
        self.about_dialog = AboutDialog()
        self.about_dialog.show()

    def show_disclaimer(self):
        """
        Opens a dialog box with some message.

        :return: None
        :rtype: NoneType

        """
        self.disclaimer_dialog = DisclaimerDialog()
        self.disclaimer_dialog.show()

    def load_sites(self) -> None:
        """
        Method to load the sites in the time_series folder in the project folder.
        It displays all possible sites for local and remote reference processing.

        :return: None
        :rtype: NoneType

        """
        # Resetting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")

        self.new_fs.hide()
        self.fft_length_dropdown.clear()
        self.parzen_radius_entry.clear()
        self.md_threshold_entry.clear()

        try:
            self.localsite_dropdown.currentIndexChanged.disconnect(self.site_dropdown_changed)
            self.remotesite_dropdown.currentIndexChanged.disconnect(self.site_dropdown_changed)
            self.sampling_frequency_dropdown.clear()
            self.localsite = None
            self.remotesite = None
        except:
            pass

        if self.project_dir is not None:
            try:
                self.allsites = phoenix_utils.load_sites(self.project_dir)
                if self.localsite_dropdown is not None:
                    self.localsite_dropdown.clear()
                self.localsite_dropdown.addItem("Select---")
                self.localsite_dropdown.addItems(self.allsites)
                # Disable the "Select" option
                self.localsite_dropdown.model().item(0).setEnabled(False)
                self.localsite_dropdown.currentIndexChanged.connect(self.site_dropdown_changed)
                if self.remotesite_dropdown is not None:
                    self.remotesite_dropdown.clear()
                self.remotesite_dropdown.addItem("Select---")
                self.remotesite_dropdown.addItems(self.allsites)
                self.remotesite_dropdown.currentIndexChanged.connect(self.site_dropdown_changed)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred while loading sites: {e}")
        else:
            QMessageBox.warning(self, "Warning", "Please create and/or open a project!")
            return

    def site_dropdown_changed(self) -> None:
        """
        Gets current text from local site and remote site dropdown menus and
        gets sampling frequencies available in both stations.
        Displays in the sampling frequency dropdown.

        :return: None
        :rtype: NoneType

        """
        self.new_fs.hide()
        self.fft_length_dropdown.clear()
        self.parzen_radius_entry.clear()
        self.md_threshold_entry.clear()

        self.localsite = self.localsite_dropdown.currentText()
        self.localsite_path = os.path.join(
            self.project_dir,
            'time_series',
            self.localsite
        )

        # Find the first valid recording folder
        recording_pattern = re.compile(
            r"^\d+_\d{4}-\d{2}-\d{2}-\d{6}$"
        )

        recording_folders = sorted(
            entry for entry in os.listdir(self.localsite_path)
            if os.path.isdir(os.path.join(self.localsite_path, entry))
            and recording_pattern.match(entry)
        )

        # Full path to the first recording folder
        self.localsite_path = os.path.join(
            self.localsite_path,
            recording_folders[0]
        )

        unique_samp_values = phoenix_utils.get_sampling_rate_list(
            recording_path=self.localsite_path)

        if self.remotesite_dropdown.currentText() == "Select---":
            self.remotesite = None
        else:
            self.remotesite = self.remotesite_dropdown.currentText()
            self.remotesite_path = os.path.join(
                self.project_dir,
                'time_series',
                self.remotesite
            )

        # Read metadata
        # Define the file path
        local_recmeta_file = os.path.join(
            self.localsite_path,
            'recmeta.json'
        )

        # Read JSON file
        with open(local_recmeta_file, 'r', encoding='utf-8') as f:
            self.recmeta_data_local = json.load(f)

        local_channel_map = self.recmeta_data_local['channel_map']['mapping']
        self.local_channel_map = {ch['tag']: ch['idx'] for ch in local_channel_map}

        # Finding overlapping measurements
        if self.remotesite is not None:
            remote_recording_folders = sorted(
                entry for entry in os.listdir(self.remotesite_path)
                if os.path.isdir(os.path.join(self.remotesite_path, entry))
                and recording_pattern.match(entry)
            )

            # Full path to the first recording folder
            self.remotesite_path = os.path.join(
                self.remotesite_path,
                remote_recording_folders[0]
            )

            remote_recmeta_file = os.path.join(
                self.remotesite_path,
                'recmeta.json'
            )

            # Read JSON file
            with open(remote_recmeta_file, 'r', encoding='utf-8') as f:
                self.recmeta_data_remote = json.load(f)

            remote_channel_map = self.recmeta_data_remote['channel_map']['mapping']
            self.remote_channel_map = {ch['tag']: ch['idx'] for ch in remote_channel_map}

            local_start_time = self.recmeta_data_local.get('start', None)
            remote_start_time = self.recmeta_data_remote.get('start', None)

            if local_start_time != remote_start_time:
                QMessageBox.warning(self, 'Warning',
                                    "Remote reference processing cannot be done as "
                                    "local and remote start time doesn't match. "
                                    "Currently, SigMT needs same start time.")
                self.remotesite = None
                self.remotesite_dropdown.setCurrentIndex(0)
                return

            unique_sampling_rates_rr = phoenix_utils.get_sampling_rate_list(
                recording_path=self.remotesite_path
            )

            unique_combined = sorted(set(unique_samp_values) | set(unique_sampling_rates_rr))

            verified_sampling_rates = []

            if unique_combined:
                for sampling_rate in unique_combined:
                    file_extension = phoenix_utils.sampling_rate_to_extension(
                        sampling_rate=int(sampling_rate)
                    )

                    time_stamp = phoenix_utils.return_overlapping_info(
                        file_extension=file_extension,
                        local_station_path=self.localsite_path,
                        remote_station_path=self.remotesite_path,
                    )

                    if time_stamp:
                        verified_sampling_rates.append(sampling_rate)

            # If no overlapping found with remote, return to single site mode
            if not verified_sampling_rates:
                QMessageBox.warning(self, 'Warning',
                                    "No overlapping time series found!")
                self.remotesite = None
                self.remotesite_dropdown.setCurrentIndex(0)
                verified_sampling_rates = unique_samp_values
        else:
            self.recmeta_data_remote = None
            self.remote_channel_map = None
            verified_sampling_rates = unique_samp_values

        if self.remotesite is None:
            verified_sampling_rates = unique_samp_values

        # Displaying in sampling frequency dropdown
        self.sampling_frequency_dropdown.clear()  # Clear existing dropdown
        self.sampling_frequency_dropdown.addItems(verified_sampling_rates)  # Add new values

    def read_ts(self) -> None:
        """
        Reads time series based from local and remote station (if requested).
        Type of reading is decided based on the sampling rate.

        :return: None
        :rtype: NoneType

        """
        # Resetting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        self.new_fs.hide()

        self.sfreq_selected = self.sampling_frequency_dropdown.currentText()

        if not self.sfreq_selected:
            QMessageBox.warning(self, 'Warning', "Please choose sampling frequency.")
            return

        if self.localsite is None:
            QMessageBox.warning(self, 'Warning', "Please choose local and/or remote site.")
            return

        if int(self.sfreq_selected) > 150:
            self.file_type = 'decimated_segmented'
            print('Decimated Segmented')
        else:
            self.file_type = 'decimated_continuous'
            print('Decimated Continuous')

        self.file_extension = phoenix_utils.sampling_rate_to_extension(
            sampling_rate=int(self.sfreq_selected)
        )

        # Reading time series data
        if self.file_type == 'decimated_continuous':
            self.time_series = phoenix_utils.read_decimated_continuous_data(
                recording_path=self.localsite_path,
                channel_map=self.local_channel_map,
                file_extension=self.file_extension
            )

            if self.remotesite:
                remote_time_series = phoenix_utils.read_decimated_continuous_data(
                    recording_path=self.remotesite_path,
                    channel_map=self.remote_channel_map,
                    file_extension=self.file_extension
                )
                remote_time_series_length = len(
                    next(iter(remote_time_series.get('run0', {}).values())))
                local_time_series_length = len(
                    next(iter(self.time_series.get('run0', {}).values())))

                min_time_series_length = min(remote_time_series_length, local_time_series_length)

                # Assuming same start time for local and remote station
                for run in remote_time_series.keys():
                    for channel in remote_time_series[run].keys():
                        remote_time_series[run][channel] = remote_time_series[run][channel][
                                                           :min_time_series_length]

                for run in self.time_series.keys():
                    for channel in self.time_series[run].keys():
                        self.time_series[run][channel] = self.time_series[run][channel][
                                                         :min_time_series_length]
                    if 'hx' in remote_time_series[run].keys():
                        self.time_series[run]['rx'] = remote_time_series[run]['hx'].copy()
                    if 'hy' in remote_time_series[run].keys():
                        self.time_series[run]['ry'] = remote_time_series[run]['hy'].copy()

        elif self.file_type == 'decimated_segmented':
            self.time_series = phoenix_utils.read_decimated_segmented_data(
                recording_path=self.localsite_path,
                channel_map=self.local_channel_map,
                file_extension=self.file_extension
            )
            if self.remotesite:
                remote_time_series = phoenix_utils.read_decimated_segmented_data(
                    recording_path=self.remotesite_path,
                    channel_map=self.remote_channel_map,
                    file_extension=self.file_extension
                )

                num_remote_segments = len(remote_time_series.keys())
                num_local_segments = len(self.time_series.keys())

                min_num_segments = min(num_remote_segments, num_local_segments)

                self.time_series = dict(list(self.time_series.items())[:min_num_segments])
                remote_time_series = dict(
                    list(remote_time_series.items())[:min_num_segments])

                for run in self.time_series.keys():
                    if 'hx' in remote_time_series[run].keys():
                        self.time_series[run]['rx'] = remote_time_series[run]['hx'].copy()
                    if 'hy' in remote_time_series[run].keys():
                        self.time_series[run]['ry'] = remote_time_series[run]['hy'].copy()

        # read data here
        self.procinfo['fs'] = int(self.sfreq_selected)

        # Read calibration data
        self.calibration_data_electric = phoenix_utils.prepare_calibration_data_electric(
            local_recmeta_data=self.recmeta_data_local,
            channel_map=self.local_channel_map
        )

        self.calibration_data_magnetic = phoenix_utils.prepare_calibration_data_magnetic(
            project_dir=self.project_dir,
            local_recmeta_data=self.recmeta_data_local,
            local_channel_map=self.local_channel_map,
            remote_recmeta_data=self.recmeta_data_remote,
            remote_channel_map=self.remote_channel_map,
        )

        num_samples = len(next(iter(next(iter(self.time_series.values())).values())))

        QMessageBox.information(
            self,
            "Information",
            "Summary\n"
            "----------\n"
            f"Sampling frequency: {self.procinfo['fs']} Hz\n"
            f"File type: {self.file_type.replace('_', ' ').title()}\n"
            f"Number of continuous samples: {num_samples}\n"
        )

        self.procinfo['nsamples_mostly'] = num_samples
        fftlength = utils.get_fftlength(self.procinfo['nsamples_mostly'])

        # Updating FFT length dropdown
        fft_values = [
            v for v in self.fft_values
            if int(v) <= self.procinfo["nsamples_mostly"]
        ]

        fft_values = [str(v) for v in fft_values if int(v) < self.procinfo["nsamples_mostly"]]

        self.fft_length_dropdown.addItems(fft_values)

        if self.file_type == 'decimated_continuous':
            self.fft_length_dropdown.setCurrentIndex(
                self.fft_length_dropdown.findText(str(fftlength)))
        elif self.file_type == 'decimated_segmented':
            self.fft_length_dropdown.setCurrentIndex(0)
        # Updating parzen window radius
        parzen_radius = utils.get_parzen(self.procinfo['fs'])
        self.parzen_radius_entry.setText(str(parzen_radius))
        # Updating MD thresholds
        self.md_threshold_entry.setText(str(1.5))

    def decimate(self) -> None:
        """
        Decimate the time series.

        :return: None
        :rtype: NoneType

        """
        # Reseting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        #
        if self.localsite is None:
            QMessageBox.warning(self, 'Warning', "Please choose local and/or remote site.")
            return
        try:
            progress_dialog = QProgressDialog("Performing Decimation...", None, 0, 1, self)
            progress_dialog.setWindowModality(Qt.WindowModal)
            progress_dialog.setWindowTitle("Please wait")
            progress_dialog.show()
            # This is important to show the progress bar and GUI not to freeze
            qapp_instance = QApplication.instance()
            qapp_instance.processEvents()
            decimation_factor = int(self.decimation_list_dropdown.currentText())

            for ts in self.time_series:
                for channel in self.time_series[ts]:
                    self.time_series[ts][channel] = signal.decimate(
                        self.time_series[ts][channel][:],
                        decimation_factor,
                        n=None,
                        ftype='iir'
                    )

            self.procinfo['fs'] = self.procinfo['fs'] / decimation_factor
            self.new_fs.setText(f"Now sampling frequency is {self.procinfo['fs']} Hz")
            self.new_fs.show()

            first_run = next(iter(self.time_series.values()))
            first_array = next(iter(first_run.values()))
            self.procinfo['nsamples_mostly'] = first_array.shape[0]

            fftlength = utils.get_fftlength(self.procinfo['nsamples_mostly'])

            # Updating FFT length dropdown
            fft_values = [
                v for v in self.fft_values
                if int(v) <= self.procinfo["nsamples_mostly"]
            ]

            # Replace existing items
            self.fft_length_dropdown.blockSignals(True)
            self.fft_length_dropdown.clear()
            self.fft_length_dropdown.addItems([str(v) for v in fft_values])

            if self.file_type == 'decimated_continuous':
                # Select fftlength if present, otherwise pick a sensible fallback
                idx = self.fft_length_dropdown.findText(str(fftlength))
                if idx >= 0:
                    self.fft_length_dropdown.setCurrentIndex(idx)
                else:
                    # fallback: last item (largest) if any
                    if self.fft_length_dropdown.count() > 0:
                        self.fft_length_dropdown.setCurrentIndex(
                            self.fft_length_dropdown.count() - 1)
            elif self.file_type == 'decimated_segmented':
                self.fft_length_dropdown.setCurrentIndex(0)

            self.fft_length_dropdown.blockSignals(False)

            # Updating parzen window radius
            parzen_radius = utils.get_parzen(self.procinfo['fs'])
            self.parzen_radius_entry.setText(str(parzen_radius))
            #
            progress_dialog.setValue(1)
            qapp_instance.processEvents()
            QMessageBox.information(self, "Done", "Decimation done!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while decimation: {e}")

    def handle_radio_button(self) -> None:
        """
        Handles the state change of the radio buttons and updates the notch filter status.

        :return: None
        :rtype: NoneType

        """
        if self.notch_radio_off.isChecked():
            self.notch_status = 'off'
        elif self.notch_radio_on.isChecked():
            self.notch_status = 'on'

    def save_parameters(self) -> None:
        """
        Saves the parameters to the self.procinfo dictionary.

        :return: None
        :rtype: NoneType

        """
        # Reseting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        #
        try:
            if self.time_series is not None:
                self.procinfo['localsite'] = self.localsite
                self.procinfo['remotesite'] = self.remotesite
                self.procinfo['notch'] = self.notch_status
                self.procinfo['processing_mode'] = self.project_setup['processing_mode']
                self.procinfo['fft_length'] = int(self.fft_length_dropdown.currentText())
                self.procinfo['parzen_radius'] = float(self.parzen_radius_entry.text())
                self.procinfo['md_thresh'] = float(self.md_threshold_entry.text())
                self.procinfo['notch_frequency'] = float(self.project_setup['notch_frequency'])
                self.procinfo['target_frequency_table_type'] = self.project_setup[
                    'target_frequency_table_type'].lower()
                self.procinfo['frequencies_per_decade'] = int(
                    self.project_setup['frequencies_per_decade'])

                self.procinfo['lat'] = self.recmeta_data_local['timing']['gps_lat']
                self.procinfo['lon'] = self.recmeta_data_local['timing']['gps_lon']
                self.procinfo['elev'] = self.recmeta_data_local['timing']['gps_alt']
                #
                self.procinfo['start_time'] = 0  # TODO
                QMessageBox.information(self, "Done", f"Parameters are saved.")
            else:
                QMessageBox.critical(self, "Error", f"Please read time series first!!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}."
                                                "Ensure time series is read!!!")

    def perform_bandavg(self) -> None:
        """
        Perform band averaging.

        :return: None
        :rtype: NoneType

        """
        # Resetting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        #
        if 'localsite' in self.procinfo and self.procinfo['localsite'] == self.localsite:

            if self.file_type == 'decimated_continuous':
                self.band_averaging_decimated_continuous()
            elif self.file_type == 'decimated_segmented':
                self.band_averaging_decimated_segmented()

            # Calculating data selection parameters
            time_dataselection = time.time()
            self.bandavg_dataset['coh_ex'] = (
                ('time_window', 'frequency'), dstools.cohex(self.bandavg_dataset))
            self.bandavg_dataset['coh_ey'] = (
                ('time_window', 'frequency'), dstools.cohey(self.bandavg_dataset))
            if not self.project_setup['processing_mode'] == "MT Only":
                self.bandavg_dataset['coh_hz'] = (
                    ('time_window', 'frequency'), dstools.cohhz(self.bandavg_dataset))
            self.bandavg_dataset['alpha_h'], self.bandavg_dataset[
                'alpha_e'] = dstools.pdvalues(
                self.bandavg_dataset)
            print(f'Time taken for data selection tool: ' + str(time.time() - time_dataselection))
            QMessageBox.information(self, "Done", f"Band averaging done!")
        else:
            QMessageBox.critical(self, "Error", f"Please read time series again!!")

    def band_averaging_decimated_continuous(self):
        datasets = []

        progress_dialog = QProgressDialog(
            "Performing band averaging...",
            None,
            0,
            len(self.time_series),
            self
        )
        progress_dialog.setWindowModality(Qt.WindowModal)
        progress_dialog.setWindowTitle("Please wait")
        progress_dialog.show()
        # This is important to show the progress bar and GUI not to freeze
        qapp_instance = QApplication.instance()
        qapp_instance.processEvents()
        #
        num = 0
        bandavg_time = time.time()

        for run in self.time_series:
            # Preparing inputs for band averaging
            # TODO: Replace this with a better strategy later
            if self.procinfo['notch'] == 'on':
                notch_filter_apply = True
            else:
                notch_filter_apply = False

            # TODO: Replace this with a better strategy later
            if self.procinfo['processing_mode'] == 'MT + Tipper':
                process_mt = True
                process_tipper = True
            elif self.procinfo['processing_mode'] == 'MT Only':
                process_mt = True
                process_tipper = False
            else:
                process_mt = None
                process_tipper = None

            # TODO: Replace this with a better strategy later
            if self.procinfo['remotesite'] is not None:
                remote_reference = True
            else:
                remote_reference = False

            # Get the bandavg object
            bandavg = BandAveraging(
                time_series=self.time_series[run].copy(),
                sampling_frequency=self.procinfo['fs'],
                fft_length=self.procinfo['fft_length'],
                parzen_window_radius=self.procinfo['parzen_radius'],
                overlap=50,
                frequencies_per_decade=self.procinfo['frequencies_per_decade'],
                remote_reference=remote_reference,
                calibrate_electric=True,
                calibrate_magnetic=True,
                calibration_data_electric=self.calibration_data_electric,
                calibration_data_magnetic=self.calibration_data_magnetic,
                instrument='phoenix',
                apply_notch_filter=notch_filter_apply,
                notch_frequency=self.procinfo['notch_frequency'],
                process_mt=process_mt,
                process_tipper=process_tipper,
            )

            datasets.append(bandavg.band_averaged_dataset)  # appends xarray dataset (for a run)
            num += 1
            progress_dialog.setValue(num)
        self.bandavg_dataset = xr.concat(datasets, dim='time_window').assign_coords(
            time_window=np.arange(len(xr.concat(datasets, dim='time_window').time_window)))
        self.dof = bandavg.dof
        self.avgf = bandavg.avgf
        self.bandavg_dataset['dof'] = xr.DataArray(
            self.dof,
            coords={'frequency': self.bandavg_dataset.coords['frequency']},
            dims='frequency'
        )
        print(f'Time taken for band averaging: ' + str(time.time() - bandavg_time))
        qapp_instance.processEvents()

    def band_averaging_decimated_segmented(self):
        datasets = []

        # decimated_segmented type usually contains thousands of
        # independent time series which make band averaging inefficient.
        # so adding an optimization technique to reduce number of runs
        optimized_time_series = phoenix_utils.optimize_time_series_dict(
            time_series_dict=self.time_series.copy(),
            fft_length=self.procinfo['fft_length'],
            overlap=50,
        )

        progress_dialog = QProgressDialog(
            "Performing band averaging...",
            None,
            0,
            len(optimized_time_series),
            self
        )
        progress_dialog.setWindowModality(Qt.WindowModal)
        progress_dialog.setWindowTitle("Please wait")
        progress_dialog.show()
        # This is important to show the progress bar and GUI not to freeze
        qapp_instance = QApplication.instance()
        qapp_instance.processEvents()
        #
        num = 0
        bandavg_time = time.time()

        for run in optimized_time_series:

            # Preparing inputs for band averaging
            # TODO: Replace this with a better strategy later
            if self.procinfo['notch'] == 'on':
                notch_filter_apply = True
            else:
                notch_filter_apply = False

            # TODO: Replace this with a better strategy later
            if self.procinfo['processing_mode'] == 'MT + Tipper':
                process_mt = True
                process_tipper = True
            elif self.procinfo['processing_mode'] == 'MT Only':
                process_mt = True
                process_tipper = False
            else:
                process_mt = None
                process_tipper = None

            # TODO: Replace this with a better strategy later
            if self.procinfo['remotesite'] is not None:
                remote_reference = True
            else:
                remote_reference = False

            # Get the bandavg object
            bandavg = BandAveraging(
                time_series=optimized_time_series[run],
                sampling_frequency=self.procinfo['fs'],
                fft_length=self.procinfo['fft_length'],
                parzen_window_radius=self.procinfo['parzen_radius'],
                overlap=50,
                frequencies_per_decade=self.procinfo['frequencies_per_decade'],
                remote_reference=remote_reference,
                calibrate_electric=True,
                calibrate_magnetic=True,
                calibration_data_electric=self.calibration_data_electric,
                calibration_data_magnetic=self.calibration_data_magnetic,
                instrument='phoenix',
                apply_notch_filter=notch_filter_apply,
                notch_frequency=self.procinfo['notch_frequency'],
                process_mt=process_mt,
                process_tipper=process_tipper,
                reshape=False,
            )

            datasets.append(bandavg.band_averaged_dataset)  # appends xarray dataset (for a run)
            num += 1
            progress_dialog.setValue(num)
        self.bandavg_dataset = xr.concat(datasets, dim='time_window').assign_coords(
            time_window=np.arange(len(xr.concat(datasets, dim='time_window').time_window)))
        self.dof = bandavg.dof
        self.avgf = bandavg.avgf
        self.bandavg_dataset['dof'] = xr.DataArray(
            self.dof,
            coords={'frequency': self.bandavg_dataset.coords['frequency']},
            dims='frequency'
        )
        print(f'Time taken for band averaging: ' + str(time.time() - bandavg_time))
        qapp_instance.processEvents()

    def plot_coh_all(self) -> None:
        """
        Plot coherency plots for all frequencies

        :return: None
        :rtype: NoneType

        """
        plots.plot_coherencies_all(self.bandavg_dataset)

    def plot_pd_all(self) -> None:
        """
        Plot polarization direction plots for all frequencies

        """
        plots.plot_pd_all(self.bandavg_dataset)

    def apply_coh_thresh(self) -> None:
        """
        Apply coherency threshold.

        :return: None
        :rtype: NoneType

        """
        if self.project_setup:
            output_channels = None
            coh_thresh = float(self.coherency_threshold_edit.text())
            min_percent = float(self.minimum_data_edit.text())
            if self.project_setup['processing_mode'] == "MT + Tipper":
                output_channels = ['ex', 'ey', 'hz']
            elif self.project_setup['processing_mode'] == "MT Only":
                output_channels = ['ex', 'ey']
            self.bandavg_dataset = pds.perform_coh_thresh(self.bandavg_dataset, coh_thresh,
                                                          min_percent,
                                                          output_channels)

            self.apply_coh_thresh_button.setText("Coherency threshold APPLIED!")

    def clear_coh_thresh(self) -> None:
        """
        Clear coherency threshold value.

        :return: None
        :rtype: NoneType

        """
        if self.bandavg_dataset:
            self.bandavg_dataset = pds.clear_coh_thresh(self.bandavg_dataset)
            self.apply_coh_thresh_button.setText("Apply coherency threshold")

    def apply_pd_thresh(self) -> None:
        """
        Apply polarization direction threshold.

        :return: None
        :rtype: NoneType

        """
        if self.bandavg_dataset:
            if not all(item in self.bandavg_dataset for item in ['ex', 'ey', 'hx', 'hy']):
                QMessageBox.warning(self, "Warning", "One of ['Ex', 'Ey', 'Hx', 'Hy'] is missing.")
            else:
                component = self.pd_combo_box.currentText()
                pd_min = float(self.pd_min_edit.text())
                pd_max = float(self.pd_max_edit.text())
                self.bandavg_dataset = pds.perform_pd_selection(self.bandavg_dataset, component,
                                                                pd_min, pd_max)
                self.apply_pd_thresh_button.setText("PD threshold APPLIED!")

    def clear_pd_thresh(self) -> None:
        """
        Clear polarization direction threshold.

        :return: None
        :rtype: NoneType

        """
        if self.bandavg_dataset:
            self.bandavg_dataset = pds.clear_pd_selection(self.bandavg_dataset)
            self.apply_pd_thresh_button.setText("Perform PD thresholding")

    def perform_robust_estimation(self) -> None:
        """
        Performs robust estimation with UI handling.

        :return: None
        :rtype: NoneType

        """
        if self.bandavg_dataset:
            progress_dialog = QProgressDialog("Performing Robust Estimation...", None, 0, 1, self)
            progress_dialog.setWindowModality(Qt.WindowModal)
            progress_dialog.setWindowTitle("Please wait")
            progress_dialog.show()
            # This is important to show the progress bar and GUI not to freeze
            qapp_instance = QApplication.instance()
            qapp_instance.processEvents()

            # Perform the robust estimation
            estimation_instance = RobustEstimation(self.procinfo, self.bandavg_dataset)
            self.estimates = estimation_instance.estimates

            progress_dialog.setValue(1)
            qapp_instance.processEvents()
            QMessageBox.information(self, "Done", "Robust estimation done!")

    def plot_impedance(self) -> None:
        """
        Plot apparent resitivity and phase.

        :return: None
        :rtype: NoneType

        """
        plots.plot_mt_app_res(self.estimates, self.procinfo)

    def plot_tipper(self) -> None:
        """
        Plot tipper.

        :return: None
        :rtype: NoneType

        """
        if self.project_setup:
            if not self.project_setup['processing_mode'] == "MT Only":
                plots.plot_tipper(self.estimates, self.procinfo)

    def plot_coherency(self) -> None:
        """
        Plot coherency.

        :return: None
        :rtype: NoneType

        """
        plots.plot_coherency(self.estimates, self.procinfo)

    def save_as_edi(self) -> None:
        """
        Save as EDI.

        :return: None
        :rtype: NoneType

        """
        edi_ops.save_edi(self.estimates, self.procinfo, self.project_setup)

    def open_edi_merger(self) -> None:
        """
        Open EDI merger.

        :return: None
        :rtype: NoneType

        """
        self.edi_merger = EDIMerger()
        self.edi_merger.show()

    def exit_all(self) -> None:
        """
        Close the application.

        This method is used to exit the application by closing the main window.
        It calls the self.close() method, which typically triggers the
        application to terminate.

        :return: None
        :rtype: NoneType

        """
        self.close()
