"""
Design of the GUI for the editing of project settings.
"""
from PyQt5.QtWidgets import QLabel, QDialog, QLineEdit, QVBoxLayout, QDialogButtonBox, QComboBox


class EditProjectSetupDialog(QDialog):
    """
    GUI to edit project setup
    """

    def __init__(self, parent=None, interface=None, project_setup=None):
        super().__init__(parent)
        self.interface = interface
        self.project_setup = project_setup
        self.setWindowTitle("Edit Project Setup")
        layout = QVBoxLayout()

        self.project_name = QLineEdit(self.project_setup['project_name'])
        self.country = QLineEdit(self.project_setup['country'])
        self.acquired_by_company = QLineEdit(self.project_setup['acquired_by'])
        #
        self.processing_mode = QComboBox()
        self.processing_mode.addItem('MT + Tipper')
        self.processing_mode.addItem('MT Only')
        # TODO: Add this in future
        # self.processing_mode.addItem('Tipper Only')
        self.processing_mode.setCurrentText(self.project_setup['processing_mode'])

        self.frequencies_per_decade = QLineEdit(self.project_setup['frequencies_per_decade'])

        self.frequency_table_type = QComboBox()
        self.frequency_table_type.addItem('Default')
        if self.interface == 'Metronix':
            self.frequency_table_type.addItem('Metronix')
        # Set combo box selection based on dict key (with fallback to 'Default')
        target_type = self.project_setup.get('target_frequency_table_type', 'Default')
        self.frequency_table_type.setCurrentText(target_type)
        self.toggle_frequencies_per_decade(text= target_type)
        self.frequency_table_type.currentTextChanged.connect(self.toggle_frequencies_per_decade)

        self.notch_frequency = QLineEdit(self.project_setup['notch_frequency'])
        self.preferred_cal_file = QComboBox()
        if self.interface == 'Metronix':
            self.preferred_cal_file.addItem('xml')
            self.preferred_cal_file.addItem('metronix_txt')
        elif self.interface == 'Phoenix':
            self.preferred_cal_file.addItem('phoenix_json')
        self.preferred_cal_file.setCurrentText(self.project_setup['preferred_cal_file'])

        layout.addWidget(QLabel("Project Name:"))
        layout.addWidget(self.project_name)

        layout.addWidget(QLabel("Country:"))
        layout.addWidget(self.country)

        layout.addWidget(QLabel("Acquired by:"))
        layout.addWidget(self.acquired_by_company)

        layout.addWidget(QLabel("Processing mode:"))
        layout.addWidget(self.processing_mode)

        layout.addWidget(QLabel("Target frequency table type:"))
        layout.addWidget(self.frequency_table_type)

        layout.addWidget(QLabel("Frequencies per Decade:"))
        layout.addWidget(self.frequencies_per_decade)

        layout.addWidget(QLabel("Notch Frequency:"))
        layout.addWidget(self.notch_frequency)

        layout.addWidget(QLabel("Preferred calibration file:"))
        layout.addWidget(self.preferred_cal_file)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.setLayout(layout)

    def get_values(self) -> dict:
        """
        Returns dictionary with values

        :return: Dictionary with values
        :rtype: dict
        """
        return {
            "project_name": self.project_name.text(),
            "country": self.country.text(),
            "acquired_by": self.acquired_by_company.text(),
            "processing_mode": self.processing_mode.currentText(),
            "target_frequency_table_type": self.frequency_table_type.currentText(),
            "frequencies_per_decade": self.frequencies_per_decade.text(),
            "notch_frequency": self.notch_frequency.text(),
            "interface": self.interface,
            "preferred_cal_file": self.preferred_cal_file.currentText()
        }

    def toggle_frequencies_per_decade(self, text):
        """
        Change frequencies per decade input based on the type of frequency table
        selected.

        :param text: The selected frequency table type from the combo box.
        :type text: str
        :return: None
        :rtype: NoneType

        """
        if text == 'Metronix':
            self.frequencies_per_decade.setText('8')
            self.frequencies_per_decade.setEnabled(False)
        else:
            self.frequencies_per_decade.setEnabled(True)
            self.frequencies_per_decade.setText('12')
