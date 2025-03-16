"""
The design of all dialog boxes appears from the metronix main window
"""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QListWidget, QVBoxLayout, QDialogButtonBox, QGroupBox, \
    QHBoxLayout, QLabel, \
    QLineEdit, QComboBox, QGridLayout


class SelectionDialog(QDialog):
    """
    There are cases where user needs to manually select the remote measurement. This class
    is the design of user interface to manually select the measurement for remote reference.
    It shows two boxes with all measurements. User can select the measurements and submit.
    """

    def __init__(self, local_meas, remote_meas, parent=None):
        """
        Constructor
        :param local_meas: Measurements in the local site
        :type local_meas: list
        :param remote_meas: Measurements in the remote site
        :type remote_meas: list
        :param parent: Parent
        :type parent: object

        """
        super().__init__(parent)

        self.setWindowTitle("Select corresponding single and remote measurements")

        # Create two list widgets within group boxes
        self.local_meas_groupbox = QGroupBox("Local Measurements")
        self.local_meas_list = QListWidget()
        self.local_meas_groupbox.setLayout(QVBoxLayout())
        self.local_meas_groupbox.layout().addWidget(self.local_meas_list)

        self.remote_meas_groupbox = QGroupBox("Remote Measurements")
        self.remote_meas_list = QListWidget()
        self.remote_meas_groupbox.setLayout(QVBoxLayout())
        self.remote_meas_groupbox.layout().addWidget(self.remote_meas_list)

        for name in local_meas:
            self.local_meas_list.addItem(name)
        for name in remote_meas:
            self.remote_meas_list.addItem(name)

        self.button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        # Set up the layout
        layout = QVBoxLayout()
        hbox = QHBoxLayout()
        hbox.addWidget(self.local_meas_groupbox)
        hbox.addWidget(self.remote_meas_groupbox)
        layout.addLayout(hbox)
        layout.addWidget(self.button_box)

        self.setLayout(layout)

    def get_selected_values(self) -> tuple:
        """
        Returns the values selected by the user.

        :return: Measurement (meas) names that user selects.
        :rtype: tuple

        """
        return tuple(
            lst.currentItem().text() if lst.currentItem() else None
            for lst in (self.local_meas_list, self.remote_meas_list)
        )


class LayoutSettingsDialog(QDialog):
    """
    It is the layout to show MT field layout settings. Presently it shows the dipole lengths and
    coil numbers only.

    """

    def __init__(self, header, parent=None):
        """
        Constructor

        :param header: Header information of all time series loaded.
        :type header: Dict
        :param parent: Parent
        :type parent: object

        """
        super().__init__(parent)
        self.header = header
        self.setWindowTitle("Layout settings")
        main_layout = QVBoxLayout(self)
        layout = QGridLayout()
        main_layout.addLayout(layout)
        list_of_coils = ['MFS06', 'MFS06e', 'MFS07e']

        if 'ex' in self.header[next(iter(self.header))]:
            self.dipole_north = QLineEdit(
                str(abs(self.header[next(iter(self.header))]['ex']['x1'][0])))
            self.dipole_south = QLineEdit(
                str(abs(self.header[next(iter(self.header))]['ex']['x2'][0])))
        if 'ey' in self.header[next(iter(self.header))]:
            self.dipole_east = QLineEdit(
                str(abs(self.header[next(iter(self.header))]['ey']['y1'][0])))
            self.dipole_west = QLineEdit(
                str(abs(self.header[next(iter(self.header))]['ey']['y1'][0])))
        #
        if 'hx' in self.header[next(iter(self.header))]:
            self.coil_hx_type = QComboBox()
            self.coil_hx_type.addItems(list_of_coils)
            initial_hx_type = str(self.header[next(iter(self.header))]['hx']['sensor']).replace(
                '\x00', '')
            if initial_hx_type in list_of_coils:
                self.coil_hx_type.setCurrentIndex(self.coil_hx_type.findText(initial_hx_type))
            self.coil_hx = QLineEdit(
                str(self.header[next(iter(self.header))]['hx']['sensor_no'][0]))
        #
        if 'hy' in self.header[next(iter(self.header))]:
            self.coil_hy_type = QComboBox()
            self.coil_hy_type.addItems(list_of_coils)
            initial_hy_type = str(self.header[next(iter(self.header))]['hy']['sensor']).replace(
                '\x00', '')
            if initial_hy_type in list_of_coils:
                self.coil_hy_type.setCurrentIndex(self.coil_hy_type.findText(initial_hy_type))
            self.coil_hy = QLineEdit(
                str(self.header[next(iter(self.header))]['hy']['sensor_no'][0]))
        #
        if 'hz' in self.header[next(iter(self.header))]:
            self.coil_hz_type = QComboBox()
            self.coil_hz_type.addItems(list_of_coils)
            initial_hz_type = str(self.header[next(iter(self.header))]['hz']['sensor']).replace(
                '\x00', '')
            if initial_hz_type in list_of_coils:
                self.coil_hz_type.setCurrentIndex(self.coil_hz_type.findText(initial_hz_type))
            self.coil_hz = QLineEdit(
                str(self.header[next(iter(self.header))]['hz']['sensor_no'][0]))
        #
        if 'rx' in self.header[next(iter(self.header))] and 'ry' in self.header[
            next(iter(self.header))]:
            #
            self.coil_rx_type = QComboBox()
            self.coil_rx_type.addItems(list_of_coils)
            initial_rx_type = str(self.header[next(iter(self.header))]['rx']['sensor']).replace(
                '\x00', '')
            if initial_rx_type in list_of_coils:
                self.coil_rx_type.setCurrentIndex(self.coil_rx_type.findText(initial_rx_type))
            self.coil_rx = QLineEdit(
                str(self.header[next(iter(self.header))]['rx']['sensor_no'][0]))
            #
            self.coil_ry_type = QComboBox()
            self.coil_ry_type.addItems(list_of_coils)
            initial_ry_type = str(self.header[next(iter(self.header))]['ry']['sensor']).replace(
                '\x00', '')
            if initial_ry_type in list_of_coils:
                self.coil_ry_type.setCurrentIndex(self.coil_ry_type.findText(initial_ry_type))
            self.coil_ry = QLineEdit(
                str(self.header[next(iter(self.header))]['ry']['sensor_no'][0]))
        else:
            self.coil_rx = QLineEdit('None')
            self.coil_ry = QLineEdit('None')

        if 'ex' in self.header[next(iter(self.header))]:
            layout.addWidget(QLabel("Dipole N (m):"), 0, 0)
            layout.addWidget(self.dipole_north, 0, 1)
            layout.addWidget(QLabel("Dipole S (m):"), 0, 2)
            layout.addWidget(self.dipole_south, 0, 3)
        if 'ey' in self.header[next(iter(self.header))]:
            layout.addWidget(QLabel("Dipole E (m):"), 1, 0)
            layout.addWidget(self.dipole_east, 1, 1)
            layout.addWidget(QLabel("Dipole W (m):"), 1, 2)
            layout.addWidget(self.dipole_west, 1, 3)
        if 'hx' in self.header[next(iter(self.header))]:
            layout.addWidget(QLabel("Coil type Hx:"), 2, 0)
            layout.addWidget(self.coil_hx_type, 2, 1)
            layout.addWidget(QLabel("Coil # Hx:"), 2, 2)
            layout.addWidget(self.coil_hx, 2, 3)
        if 'hy' in self.header[next(iter(self.header))]:
            layout.addWidget(QLabel("Coil type Hy:"), 3, 0)
            layout.addWidget(self.coil_hy_type, 3, 1)
            layout.addWidget(QLabel("Coil # Hy:"), 3, 2)
            layout.addWidget(self.coil_hy, 3, 3)
        if 'hz' in self.header[next(iter(self.header))]:
            layout.addWidget(QLabel("Coil type Hz:"), 4, 0)
            layout.addWidget(self.coil_hz_type, 4, 1)
            layout.addWidget(QLabel("Coil # Hz:"), 4, 2)
            layout.addWidget(self.coil_hz, 4, 3)

        if 'rx' in self.header[next(iter(self.header))] and 'ry' in self.header[
            next(iter(self.header))]:
            layout.addWidget(QLabel("Coil type Rx:"), 5, 0)
            layout.addWidget(self.coil_rx_type, 5, 1)

            layout.addWidget(QLabel("Coil # Rx:"), 5, 2)
            layout.addWidget(self.coil_rx, 5, 3)
            #
            layout.addWidget(QLabel("Coil type Ry:"), 6, 0)
            layout.addWidget(self.coil_ry_type, 6, 1)

            layout.addWidget(QLabel("Coil # Ry:"), 6, 2)
            layout.addWidget(self.coil_ry, 6, 3)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept_vals)
        buttons.rejected.connect(self.reject)
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(buttons)

        main_layout.addLayout(buttons_layout)
        main_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

    def accept_vals(self) -> None:
        """
        Update the header with the values from the input fields

        """
        for key in self.header:
            self.header[key]['ex']['x1'][0] = float(self.dipole_north.text())
            self.header[key]['ex']['x2'][0] = float(self.dipole_south.text())
            self.header[key]['ey']['y1'][0] = float(self.dipole_east.text())
            self.header[key]['ey']['y2'][0] = float(self.dipole_west.text())
            self.header[key]['hx']['sensor'] = self.coil_hx_type.currentText()
            self.header[key]['hx']['sensor_no'][0] = int(self.coil_hx.text())
            self.header[key]['hy']['sensor'] = self.coil_hy_type.currentText()
            self.header[key]['hy']['sensor_no'][0] = int(self.coil_hy.text())
            self.header[key]['hz']['sensor'] = self.coil_hz_type.currentText()
            self.header[key]['hz']['sensor_no'][0] = int(self.coil_hz.text())
            if 'rx' in self.header[key] and 'ry' in self.header[key]:
                self.header[key]['rx']['sensor'] = self.coil_rx_type.currentText()
                self.header[key]['rx']['sensor_no'][0] = int(self.coil_rx.text())
                self.header[key]['ry']['sensor'] = self.coil_ry_type.currentText()
                self.header[key]['ry']['sensor_no'][0] = int(self.coil_ry.text())
        super().accept()
