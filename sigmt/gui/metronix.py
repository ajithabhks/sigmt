"""
First landing window for the Metronix specific operations.

"""

import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (QMainWindow, QAction, QFileDialog, QMessageBox,
                             QWidget, QVBoxLayout, QComboBox, QLabel,
                             QPushButton, QHBoxLayout, QRadioButton, QGroupBox,
                             QGridLayout, QLineEdit, QSizePolicy, QDialog,
                             QApplication, QProgressDialog)
from scipy import signal

import sigmt.utils.metronix.cal_from_metronix_txt
from sigmt.core import data_selection_tools as dstools
from sigmt.core import perform_data_selection as pds
from sigmt.core import plots
from sigmt.core.band_averaging import BandAveraging
from sigmt.core.robust_estimation import RobustEstimation
from sigmt.gui.about_dialog import AboutDialog
from sigmt.gui.edi_merger import EDIMerger
from sigmt.gui.metronix_dialogs import LayoutSettingsDialog
from sigmt.gui.metronix_dialogs import SelectionDialog
from sigmt.gui.project_related.create_project import ProjectSetupDialog
from sigmt.gui.project_related.edit_project import EditProjectSetupDialog
from sigmt.utils import utils
from sigmt.utils.edi import edi_ops
from sigmt.utils.metronix import metronix_utils

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
        self.add_remote_manual_button = None
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
        self.overlapping_meas = None
        self.sfreq_selected = None
        self.processing_df = None
        self.header = None
        self.xml_caldata = None
        self.menubar = None
        self.edi_merger = None
        self.plot_edi_button = None
        self.estimates = None
        self.perform_robust_estimation_button = None
        self.bandavg_dataset = None
        self.procinfo = {}
        self.h5file = None
        self.project_dir = None
        self.project_setup = None
        self.localsite = None
        self.localsite_meas = None
        self.localsite_dropdown = None
        self.remotesite_dropdown = None
        self.sampling_frequency_dropdown = None
        self.verify_layout_button = None
        self.fft_values = ['262144', '131072', '65536', '32768', '16384',
                           '8192', '4096', '2048', '1024', '512', '256']
        #
        self.remotesite = None
        self.remotesite_meas = None
        self.remotesite_manual = None
        #
        self.processing_route = None
        #
        self.interface = 'Metronix'
        self.init_ui()

    def init_ui(self) -> None:
        """
        Initial display user interface

        :return: None
        :rtype: NoneType

        """
        self.setWindowTitle('[No Project Opened] SigMT | A Tool for '
                            'Magnetotelluric Data Processing (Metronix)')
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
        remote_manual_widget = QWidget()
        remote_manual_layout = QHBoxLayout()
        remote_manual_layout.setContentsMargins(0, 0, 0, 0)
        remote_manual_widget.setLayout(remote_manual_layout)
        self.add_remote_manual_button = QPushButton("Set Remote Measurement Manually")
        self.add_remote_manual_button.clicked.connect(self.open_manual_selection_dialog)
        remote_manual_layout.addWidget(self.add_remote_manual_button)
        remote_manual_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        section1_layout.addWidget(remote_manual_widget)
        #
        section1_grid_widget = QWidget()
        section1_grid_layout = QGridLayout()
        section1_grid_layout.setContentsMargins(0, 0, 0, 0)
        section1_grid_widget.setLayout(section1_grid_layout)
        section1_grid_layout.addWidget(QLabel("Select Sampling Frequency"), 0, 0)
        self.sampling_frequency_dropdown = QComboBox()
        section1_grid_layout.addWidget(self.sampling_frequency_dropdown, 0, 1, 1, 40)
        self.read_time_series_button = QPushButton("Read time series")
        self.read_time_series_button.clicked.connect(self.read_ts)
        section1_grid_layout.addWidget(self.read_time_series_button, 1, 0)
        self.verify_layout_button = QPushButton("Verify/Edit layout settings")
        self.verify_layout_button.clicked.connect(self.open_layout_settings)
        self.verify_layout_button.hide()
        section1_grid_layout.addWidget(self.verify_layout_button, 1, 1)
        section1_grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        #
        section1_layout.addWidget(section1_grid_widget)
        section1_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

        # Subsection 2
        section2 = QGroupBox("Decimation")
        section2_layout = QGridLayout()
        section2_layout.addWidget(QLabel("Select Decimation from the list:"), 0, 0)
        self.decimation_list_dropdown = QComboBox()
        self.decimation_list_dropdown.addItems(['4', '8'])
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
                                    'A Tool for Magnetotelluric Data Processing (Metronix)')

    def open_project(self) -> None:
        """
        To open an existing project.
        
        :return: None
        :rtype: NoneType

        """
        if self.verify_layout_button is not None:
            self.verify_layout_button.hide()
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
                self.setWindowTitle(self.project_setup[
                                        'project_name'] +
                                    ' - SigMT | '
                                    'A Tool for Magnetotelluric Data Processing (Metronix)')
            else:
                QMessageBox.warning(self, "Warning", "Not a Metronix Project")
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
                self.setWindowTitle(self.project_setup[
                                        'project_name'] +
                                    ' - SigMT | '
                                    'A Tool for Magnetotelluric Data Processing (Metronix)')
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
        #
        if self.verify_layout_button is not None:
            self.verify_layout_button.hide()
        self.new_fs.hide()
        self.fft_length_dropdown.clear()
        self.parzen_radius_entry.clear()
        self.md_threshold_entry.clear()

        try:
            self.localsite_dropdown.currentIndexChanged.disconnect(self.site_dropdown_changed)
            self.remotesite_dropdown.currentIndexChanged.disconnect(self.site_dropdown_changed)
            self.sampling_frequency_dropdown.clear()
            self.localsite = None
            self.localsite_meas = None
            self.remotesite = None
            self.remotesite_meas = None
            self.remotesite_manual = None
        except:
            pass

        if self.project_dir is not None:
            try:
                self.allsites = metronix_utils.loadsites(self.project_dir)
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
                # Disable the "Select" option
                # self.remotesite_dropdown.model().item(0).setEnabled(False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred while loading sites: {e}")
        else:
            QMessageBox.warning(self, "Warning", "Please create and/or open a project!")
            return

    def site_dropdown_changed(self) -> None:
        """
        Gets current text from local site and remote site dropdown menus and finds matching
        meas folders in both. Gets sampling frequencies available in overlapping measurements.
        Displays in the sampling frequency dropdown.
        Two things are checked:
        1. If folder names are same - First and must condition
        2. Folders with same names have same type of measurements? (fs and chopper status)
        If two above conditions are met, the folders are considered as overlapping measurements.
        Added: If folder names are not matching, user will be asked to select remote meas manually.

        :TODO Need to improve strategy later. Works for now!

        :return: None
        :rtype: NoneType

        """
        if self.verify_layout_button is not None:
            self.verify_layout_button.hide()
        self.new_fs.hide()
        self.fft_length_dropdown.clear()
        self.parzen_radius_entry.clear()
        self.md_threshold_entry.clear()
        #
        self.remotesite_manual = None
        self.localsite = self.localsite_dropdown.currentText()
        if self.remotesite_dropdown.currentText() == "Select---":
            self.remotesite = None
            self.remotesite_meas = None
        else:
            self.remotesite = self.remotesite_dropdown.currentText()
        # List meas folders in local site
        self.localsite_path = os.path.join(self.project_dir, 'time_series', self.localsite)
        self.localsite_meas = metronix_utils.list_meas_folders(self.localsite_path)
        # List meas folders in remote site
        if self.remotesite is not None:
            self.remotesite_path = os.path.join(self.project_dir, 'time_series', self.remotesite)
            self.remotesite_meas = metronix_utils.list_meas_folders(self.remotesite_path)
        # Finding overlapping measurements
        if self.remotesite is not None:
            overlapping_meas = sorted(
                list(set(self.localsite_meas).intersection(set(self.remotesite_meas))))
            # If no overlapping found with remote, return to single site mode
            if not overlapping_meas:
                QMessageBox.warning(self, 'Warning',
                                    "No overlapping measurements found! You need to "
                                    "set remote measurement manually "
                                    "by clicking 'Set Remote Measurement Manually' button below. \n\n"
                                    "Else, please proceed with single site processing.")
                # self.remotesite_dropdown.blockSignals(True)
                # self.remotesite_dropdown.setCurrentIndex(0)
                # self.remotesite_dropdown.blockSignals(False)
                self.remotesite_manual = self.remotesite
                self.remotesite = None
                overlapping_meas = self.localsite_meas
        elif self.remotesite is None:
            overlapping_meas = self.localsite_meas
        # Operations on local site
        local_meas_paths = [os.path.join(self.project_dir, 'time_series', self.localsite, meas) for
                            meas in
                            overlapping_meas]
        [sampfreq, chopper_value] = metronix_utils.get_sampling_frequency_from_xml(local_meas_paths)
        localsite_meas_overlap = set(zip(overlapping_meas, sampfreq, chopper_value))
        # Operations on remote site
        if self.remotesite is not None:
            remote_meas_paths = [
                os.path.join(self.project_dir, 'time_series', self.remotesite, meas) for meas in
                overlapping_meas]
            [sampfreq, chopper_value] = metronix_utils.get_sampling_frequency_from_xml(
                remote_meas_paths)
            remotesite_meas_overlap = set(zip(overlapping_meas, sampfreq, chopper_value))
        if self.remotesite is not None:
            self.overlapping_meas = list(
                localsite_meas_overlap.intersection(remotesite_meas_overlap))
        elif self.remotesite is None:
            self.overlapping_meas = list(localsite_meas_overlap)
        # Creating processing route
        # It is a dataframe containing measurement details
        if self.remotesite is not None:
            processing_route = {
                'local': [meas[0] for meas in self.overlapping_meas],
                'remote': [meas[0] for meas in self.overlapping_meas],
                'sampling_frequency': [meas[1] for meas in self.overlapping_meas],
                'chopper_status': [meas[2] for meas in self.overlapping_meas],
                'sampling_chopper': [str(meas[1]) + ' Hz - ' + meas[2] for meas in
                                     self.overlapping_meas]
            }
        elif self.remotesite is None:
            processing_route = {
                'local': [meas[0] for meas in self.overlapping_meas],
                'remote': [None for meas in self.overlapping_meas],
                'sampling_frequency': [meas[1] for meas in self.overlapping_meas],
                'chopper_status': [meas[2] for meas in self.overlapping_meas],
                'sampling_chopper': [str(meas[1]) + ' Hz - ' + meas[2] for meas in
                                     self.overlapping_meas]
            }
        if self.remotesite_manual is not None:
            processing_route = {
                'local': [],
                'remote': [],
                'sampling_frequency': [],
                'chopper_status': [],
                'sampling_chopper': []
            }
        self.processing_route = pd.DataFrame(processing_route)
        unique_samp_values = self.processing_route['sampling_chopper'].unique().tolist()
        # Displaying in sampling frequency dropdown
        self.sampling_frequency_dropdown.clear()  # Clear existing dropdown
        self.sampling_frequency_dropdown.addItems(unique_samp_values)  # Add new values

    def open_manual_selection_dialog(self) -> None:
        """
        Opens a dialog window for manual selection of remote measurements.
        This method creates a dialog window for manually selecting remote measurements.

        :return: None
        :rtype: NoneType

        """
        try:
            dialog = SelectionDialog(local_meas=self.localsite_meas,
                                     remote_meas=self.remotesite_meas, parent=self)
            if dialog.exec_() == QDialog.Accepted:
                [local, remote] = dialog.get_selected_values()
                if self.remotesite is None:
                    self.remotesite = self.remotesite_manual
                local_meas_path = [
                    os.path.join(self.project_dir, 'time_series', self.localsite, local)]
                [l_sampfreq, l_chopper_value] = metronix_utils.get_sampling_frequency_from_xml(
                    local_meas_path)
                remote_meas_path = [
                    os.path.join(self.project_dir, 'time_series', self.remotesite, remote)]
                [r_sampfreq, r_chopper_value] = metronix_utils.get_sampling_frequency_from_xml(
                    remote_meas_path)
                if l_sampfreq != r_sampfreq:
                    QMessageBox.critical(self, "Error",
                                         f"Sampling frequencies are not matching. Cannot proceed.")
                    self.remotesite = None
                    self.remotesite_dropdown.setCurrentIndex(0)
                    return
                if l_chopper_value != r_chopper_value:
                    QMessageBox.critical(self, "Error",
                                         f"Chopper status not matching. Cannot proceed. Choose another match.")
                    self.remotesite = None
                    self.remotesite_dropdown.setCurrentIndex(0)
                    return
                new_row = {
                    'local': local,
                    'remote': remote,
                    'sampling_frequency': r_sampfreq[0],
                    'chopper_status': r_chopper_value[0],
                    'sampling_chopper': str(r_sampfreq[0]) + ' Hz - ' + r_chopper_value[0]
                }
                # Append the new row to the DataFrame
                self.processing_route = self.processing_route._append(new_row, ignore_index=True)
                unique_samp_values = self.processing_route['sampling_chopper'].unique().tolist()
                self.sampling_frequency_dropdown.clear()
                self.sampling_frequency_dropdown.addItems(unique_samp_values)
            self.processing_route.drop_duplicates(inplace=True)
        except:
            QMessageBox.critical(self, "Error",
                                 f"Select local and remote site first! It works only if remote site is selected.")

    # noinspection PyTypeChecker
    def read_ts(self) -> None:
        """
        Reads time series based on the self.processing_route.
        processing_route is a DataFrame containing ['local', 'remote', 'sampling_frequency', 'chopper_status',
        'sampling_chopper'] for all measurements.
        processing_df contains the values for selected sampling frequency & chopper status combo.

        :return: None
        :rtype: NoneType

        """
        # Reseting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        #
        self.new_fs.hide()
        if self.processing_route is not None:
            self.sfreq_selected = self.sampling_frequency_dropdown.currentText()
            self.processing_df = self.processing_route[
                self.processing_route['sampling_chopper'] == self.sfreq_selected]
        else:
            QMessageBox.warning(self, 'Warning', "Please choose sampling frequency.")
            return
        if self.localsite is None:
            QMessageBox.warning(self, 'Warning', "Please choose local and/or remote site.")
            return

        # Create a progress dialog
        progress_dialog = QProgressDialog("Reading time series...", None, 0,
                                          len(self.processing_df.index), self)
        progress_dialog.setWindowModality(Qt.WindowModal)
        progress_dialog.setWindowTitle("Please wait")
        progress_dialog.show()

        # This is important to show the progress bar and GUI not to freeze
        qapp_instance = QApplication.instance()
        qapp_instance.processEvents()

        # If single site processing or remote site processing?
        if self.remotesite is None:
            local_paths = [os.path.join(self.project_dir, 'time_series', self.localsite, meas) for
                           meas in
                           self.processing_df['local']]
            num = 0
            self.header = {}
            self.xml_caldata = {}
            self.h5file = os.path.join(self.project_dir, 'db.h5')
            if os.path.exists(self.h5file):
                os.remove(self.h5file)
            # Creates a h5 file in project directory in write mode
            with h5py.File(self.h5file, 'w') as f:
                for local_path in local_paths:
                    self.header[f'ts_{num}'], ts_dict = metronix_utils.read_ts(local_path,
                                                                               self.project_setup)
                    self.xml_caldata[f'ts_{num}'] = metronix_utils.read_calibration_from_xml(local_path)
                    ts = f.create_group(f'ts_{num}')
                    for key in ts_dict.keys():
                        ts.create_dataset(key, data=ts_dict[key].values)
                        self.header[f'ts_{num}'][key]['time_coord'] = ts_dict[key].time
                    del ts_dict
                    num += 1
                    progress_dialog.setValue(num)
        elif self.remotesite is not None:
            local_paths = [os.path.join(self.project_dir, 'time_series', self.localsite, meas) for
                           meas in
                           self.processing_df['local']]
            remote_paths = [os.path.join(self.project_dir, 'time_series', self.remotesite, meas) for
                            meas in
                            self.processing_df['remote']]
            num = 0
            self.header = {}
            self.xml_caldata = {}
            self.h5file = os.path.join(self.project_dir, 'db.h5')
            if os.path.exists(self.h5file):
                os.remove(self.h5file)
            # Creates a h5 file in project directory in write mode
            with h5py.File(self.h5file, 'w') as f:
                for local_path, remote_path in zip(local_paths, remote_paths):
                    self.header[f'ts_{num}'], ts_dict = metronix_utils.read_ts(local_path,
                                                                               self.project_setup)
                    self.xml_caldata[f'ts_{num}'] = metronix_utils.read_calibration_from_xml(local_path)

                    header_r, ts_r = metronix_utils.read_ts(remote_path, self.project_setup)

                    # Take only 'hx' and 'hy' component of remote dataset
                    header_r = {k: v for k, v in header_r.items() if k in ["hx", "hy"]}
                    ts_r = {k: v for k, v in ts_r.items() if k in ["hx", "hy"]}

                    xml_caldata_r = metronix_utils.read_calibration_from_xml(remote_path)
                    self.xml_caldata[f'ts_{num}'].update(xml_caldata_r)

                    # Some checks on data.
                    # Check if local stations time coordinates are not aligned
                    time_coords = [v.coords['time'] for v in ts_dict.values()]
                    if not all(time_coords[0].equals(tc) for tc in time_coords[1:]):
                        raise ValueError("Time coordinates are not aligned!")
                    # Check if remote stations time coordinates are not aligned
                    time_coords = [v.coords['time'] for v in ts_r.values()]
                    if not all(time_coords[0].equals(tc) for tc in time_coords[1:]):
                        raise ValueError("Time coordinates are not aligned!")

                    # Write local data to database
                    ts = f.create_group(f'ts_{num}')
                    for key in ts_dict.keys():
                        ts.create_dataset(key, data=ts_dict[key].values)

                    # Finding common time between local and remote station
                    common_time = xr.align(ts_dict[list(ts_dict.keys())[0]], ts_r[list(ts_r.keys())[0]], join="inner")[
                        0].time

                    # Write remote data to database
                    # TODO: In future, ask user if they need to user common time or use local station data
                    # TODO: where there is no overlap.
                    if len(common_time) == 0:
                        # Continue as local station processing
                        print("No overlapping time series.")
                        print("Continuing as local station processing.")
                        # Rx
                        self.header[f'ts_{num}']['rx'] = self.header[f'ts_{num}']['hx']
                        ts.create_dataset('rx', data=ts_dict['hx'].values)
                        # Ry
                        self.header[f'ts_{num}']['ry'] = self.header[f'ts_{num}']['hy']
                        ts.create_dataset('ry', data=ts_dict['hy'].values)
                    else:
                        # Else, use the overlapping time series and fill remote channels with local
                        # station data for the rest of the time
                        if 'hx' in header_r:
                            rx = ts_dict['hx'].copy()
                            rx.loc[dict(time=common_time)] = ts_r['hx'].sel(time=common_time)
                            ts.create_dataset('rx', data=rx.values)
                            header_r['hx']['nsamples'] = len(ts_r['hx'])
                            self.header[f'ts_{num}']['rx'] = header_r['hx']

                        if 'hy' in header_r:
                            ry = ts_dict['hy'].copy()
                            ry.loc[dict(time=common_time)] = ts_r['hy'].sel(time=common_time)
                            ts.create_dataset('ry', data=ry.values)
                            header_r['hy']['nsamples'] = len(ts_r['hy'])
                            self.header[f'ts_{num}']['ry'] = header_r['hy']
                    del ts_dict
                    num += 1
                    progress_dialog.setValue(num)
        qapp_instance.processEvents()
        progress_dialog.close()
        self.procinfo['fs'] = self.processing_df['sampling_frequency'].iloc[0]
        QMessageBox.information(self, "Information", "Time series reading completed!\n\n"
                                                     "Summary\n"
                                                     "----------\n"
                                                     f"No. of measurements loaded: {len(self.processing_df)}\n"
                                                     f"Sampling frequency: {self.procinfo['fs']} Hz")
        self.procinfo['nsamples_mostly'] = utils.get_nsamples(self.header)
        fftlength = utils.get_fftlength(self.procinfo['nsamples_mostly'])
        # Updating FFT length dropdown
        self.fft_length_dropdown.addItems(self.fft_values)
        self.fft_length_dropdown.setCurrentIndex(self.fft_length_dropdown.findText(str(fftlength)))
        # Updating parzen window radius
        parzen_radius = utils.get_parzen(self.procinfo['fs'])
        self.parzen_radius_entry.setText(str(parzen_radius))
        self.verify_layout_button.show()
        # Updating MD thresholds
        self.md_threshold_entry.setText(str(1.5))

    def open_layout_settings(self) -> None:
        """
        Opens layout settings dialog

        :return: None
        :rtype: NoneType

        """
        try:
            dialog = LayoutSettingsDialog(header=self.header, parent=self)
            if dialog.exec_() == QDialog.Accepted:
                self.header = dialog.header
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while layout settings: {e}")

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
            with h5py.File(self.h5file, 'r+') as f:
                for ts in f.keys():
                    for channel in f[ts]:
                        decimated_data = signal.decimate(f[ts][channel][:], decimation_factor,
                                                         n=None, ftype='iir')
                        del f[ts][channel]
                        ts_group = f[ts]
                        ts_group.create_dataset(channel, data=decimated_data)
                        self.header[ts][channel]['sfreq'] = [
                            self.header[ts][channel]['sfreq'][0] / decimation_factor]
                        self.header[ts][channel]['nsamples'] = len(decimated_data)
                        # self.header[ts][channel]['time_coord'] = self.header[ts][channel][
                        #                                              'time_coord'][
                        #                                          ::decimation_factor]
            self.procinfo['fs'] = self.procinfo['fs'] / decimation_factor
            self.new_fs.setText(f"Now sampling frequency is {self.procinfo['fs']} Hz")
            self.new_fs.show()
            self.procinfo['nsamples_mostly'] = utils.get_nsamples(self.header)
            fftlength = utils.get_fftlength(self.procinfo['nsamples_mostly'])
            # Updating FFT length dropdown
            self.fft_length_dropdown.setCurrentIndex(
                self.fft_length_dropdown.findText(str(fftlength)))
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
            if self.header is not None:
                self.procinfo['localsite'] = self.localsite
                self.procinfo['remotesite'] = self.remotesite
                self.procinfo['notch'] = self.notch_status
                self.procinfo['processing_mode'] = self.project_setup['processing_mode']
                self.procinfo['fft_length'] = int(self.fft_length_dropdown.currentText())
                self.procinfo['parzen_radius'] = float(self.parzen_radius_entry.text())
                self.procinfo['md_thresh'] = float(self.md_threshold_entry.text())
                self.procinfo['notch_frequency'] = float(self.project_setup['notch_frequency'])
                self.procinfo['preferred_cal_file'] = self.project_setup['preferred_cal_file']
                self.procinfo['target_frequency_table_type'] = self.project_setup['target_frequency_table_type'].lower()
                self.procinfo['frequencies_per_decade'] = int(
                    self.project_setup['frequencies_per_decade'])
                first_header = next(iter(next(iter(self.header.values())).values()))
                #
                self.procinfo['lat'] = first_header['lat'][0] / 1000 / 60 / 60
                self.procinfo['lon'] = first_header['lon'][0] / 1000 / 60 / 60
                self.procinfo['elev'] = first_header['elev'][0] / 100
                #
                self.procinfo['start_time'] = first_header['start'][0]
                QMessageBox.information(self, "Done", f"Parameters are saved.")
            else:
                QMessageBox.critical(self, "Error", f"Please read time series first!!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}."
                                                "Ensure time series is read!!!")

    def perform_bandavg(self) -> None:
        """
        Set up the band averaging.

        :return: None
        :rtype: NoneType

        """
        # Resetting some buttons
        self.apply_coh_thresh_button.setText("Apply coherency threshold")
        self.apply_pd_thresh_button.setText("Perform PD thresholding")
        #
        if 'localsite' in self.procinfo and self.procinfo['localsite'] == self.localsite:
            datasets = []
            #
            progress_dialog = QProgressDialog("Performing band averaging...", None, 0,
                                              len(self.header), self)
            progress_dialog.setWindowModality(Qt.WindowModal)
            progress_dialog.setWindowTitle("Please wait")
            progress_dialog.show()
            # This is important to show the progress bar and GUI not to freeze
            qapp_instance = QApplication.instance()
            qapp_instance.processEvents()
            #
            num = 0
            bandavg_time = time.time()
            for ts in self.header:

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

                # TODO: Replace this with a better strategy later
                calibration_data_electric = {}
                if 'ex' in self.header[ts]:
                    calibration_data_electric['ex'] = {}
                    calibration_data_electric['ex']['x1'] = self.header[ts]['ex']['x1'][0]
                    calibration_data_electric['ex']['x2'] = self.header[ts]['ex']['x2'][0]

                if 'ey' in self.header[ts]:
                    calibration_data_electric['ey'] = {}
                    calibration_data_electric['ey']['y1'] = self.header[ts]['ey']['y1'][0]
                    calibration_data_electric['ey']['y2'] = self.header[ts]['ey']['y2'][0]

                # TODO: Replace this with a better strategy later
                calibration_data_magnetic = {'instrument': 'metronix'}
                possible_magnetic_channels = ['hx', 'hy', 'hz', 'rx', 'ry']
                available_magnetic_channels = [element for element in possible_magnetic_channels if
                                               element in list(self.header[ts].keys())]
                for magnetic_channel in available_magnetic_channels:
                    calibration_data_magnetic[magnetic_channel] = {}
                    calibration_data_magnetic[magnetic_channel]['sensor_type'] = self.header[ts][magnetic_channel]['sensor']
                    calibration_data_magnetic[magnetic_channel]['sensor_serial_number'] = \
                        self.header[ts][magnetic_channel]['sensor_no'][0]
                    calibration_data_magnetic[magnetic_channel]['chopper_status'] = \
                        self.header[ts][magnetic_channel]['bychopper'][0]

                    if calibration_data_magnetic[magnetic_channel]['chopper_status'] == 0:
                        chopper_status = 'chopper_off'
                    elif calibration_data_magnetic[magnetic_channel]['chopper_status'] == 1:
                        chopper_status = 'chopper_on'
                    else:
                        chopper_status = None

                    print('\n')
                    print('====================================================')
                    print(f'Working on calibration data for {magnetic_channel}')
                    print(f"Coil serial number: {calibration_data_magnetic[magnetic_channel]['sensor_serial_number']}")
                    print('\n')

                    if str(calibration_data_magnetic[magnetic_channel]['sensor_serial_number']) in self.xml_caldata[
                        ts].keys():
                        cal_data_xml = self.xml_caldata[ts][
                            str(calibration_data_magnetic[magnetic_channel]['sensor_serial_number'])]
                    else:
                        cal_data_xml = None

                    # Sometimes null byte can be there in the sensor type string.
                    # Removing that from string if exists
                    if "\x00" in calibration_data_magnetic[magnetic_channel]["sensor_type"]:
                        calibration_data_magnetic[magnetic_channel]["sensor_type"] = \
                            calibration_data_magnetic[magnetic_channel]["sensor_type"].replace("\x00", "")

                    metronix_txt_filename = (calibration_data_magnetic[magnetic_channel]["sensor_type"].lower()
                                             + "_"
                                             + str(calibration_data_magnetic[magnetic_channel]["sensor_serial_number"])
                                             )
                    cal_txt_path = Path(self.project_dir) / "calibration_files" / f"{metronix_txt_filename}.txt"

                    if cal_txt_path.exists():
                        print('Metronix txt calibration file found.')
                        cal_data_txt = sigmt.utils.metronix.cal_from_metronix_txt.read_calibration_metronix_txt(
                            filepath=cal_txt_path)
                    else:
                        cal_data_txt = None
                        print(f'Metronix txt calibration file not found at {cal_txt_path}')

                    if self.project_setup['preferred_cal_file'] == 'xml':
                        print('Preferred calibration file selected: xml')
                        if (cal_data_xml is not None) and (cal_data_xml[chopper_status].size != 0):
                            print('Using calibration data from XML file.')
                            calibration_data_magnetic[magnetic_channel]['calibration_data'] = cal_data_xml
                        else:
                            print('No calibration data is available from XML file.')
                            print('Trying metronix txt.')
                            if (cal_data_txt is None) or (cal_data_txt[chopper_status].size == 0):
                                print('No calibration data is available from Metronix txt.')
                                calibration_data_magnetic[magnetic_channel]['calibration_data'] = None
                            else:
                                print('Calibration data found from Metronix txt. Using it.')
                                calibration_data_magnetic[magnetic_channel]['calibration_data'] = cal_data_txt
                    elif self.project_setup['preferred_cal_file'] == 'metronix_txt':
                        print('Preferred calibration file selected: metronix_txt')
                        if (cal_data_txt is not None) and (cal_data_txt[chopper_status].size != 0):
                            print('Using calibration data from metronix_txt file.')
                            calibration_data_magnetic[magnetic_channel]['calibration_data'] = cal_data_txt
                        else:
                            print('No calibration data is available from metronix_txt file.')
                            print('Trying XML.')
                            if cal_data_xml[chopper_status].size == 0:
                                print('No calibration data is available from XML.')
                                calibration_data_magnetic[magnetic_channel]['calibration_data'] = None
                            else:
                                print('Calibration data found from XML. Using it.')
                                calibration_data_magnetic[magnetic_channel]['calibration_data'] = cal_data_xml
                    else:
                        calibration_data_magnetic[magnetic_channel]['calibration_data'] = None
                    print('====================================================')
                    print('\n')

                # Get the bandavg object
                bandavg = BandAveraging(time_series=metronix_utils.prepare_ts_from_h5(self.h5file, ts),
                                        sampling_frequency=self.procinfo['fs'], overlap=50,
                                        calibrate_electric=True, calibrate_magnetic=True,
                                        calibration_data_electric=calibration_data_electric,
                                        calibration_data_magnetic=calibration_data_magnetic,
                                        fft_length=self.procinfo['fft_length'],
                                        parzen_window_radius=self.procinfo['parzen_radius'],
                                        target_frequency_table_type=self.procinfo['target_frequency_table_type'],
                                        frequencies_per_decade=self.procinfo['frequencies_per_decade'],
                                        apply_notch_filter=notch_filter_apply,
                                        notch_frequency=self.procinfo['notch_frequency'],
                                        process_mt=process_mt, process_tipper=process_tipper,
                                        remote_reference=remote_reference
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
            qapp_instance.processEvents()
            QMessageBox.information(self, "Done", f"Band averaging done!")
        else:
            QMessageBox.critical(self, "Error", f"Please read time series again!!")

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
