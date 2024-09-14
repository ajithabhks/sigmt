"""
Design of the GUI for the editing of project settings.
"""
from PyQt5.QtWidgets import QLabel, QDialog, QLineEdit, QVBoxLayout, QDialogButtonBox, QComboBox


class EditProjectSetupDialog(QDialog):
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
        self.processing_mode.addItem('Tipper Only')
        self.processing_mode.setCurrentText(self.project_setup['processing_mode'])
        #
        self.frequencies_per_decade = QLineEdit(self.project_setup['frequencies_per_decade'])
        self.notch_frequency = QLineEdit(self.project_setup['notch_frequency'])
        self.preferred_cal_file = QComboBox()
        self.preferred_cal_file.addItem('xml')
        self.preferred_cal_file.addItem('metronix txt')
        self.preferred_cal_file.setCurrentText(self.project_setup['preferred_cal_file'])

        layout.addWidget(QLabel("Project Name:"))
        layout.addWidget(self.project_name)

        layout.addWidget(QLabel("Country:"))
        layout.addWidget(self.country)

        layout.addWidget(QLabel("Acquired by:"))
        layout.addWidget(self.acquired_by_company)

        layout.addWidget(QLabel("Processing mode:"))
        layout.addWidget(self.processing_mode)

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

    def get_values(self):
        return {
            "project_name": self.project_name.text(),
            "country": self.country.text(),
            "acquired_by": self.acquired_by_company.text(),
            "processing_mode": self.processing_mode.currentText(),
            "frequencies_per_decade": self.frequencies_per_decade.text(),
            "notch_frequency": self.notch_frequency.text(),
            "interface": self.interface,
            "preferred_cal_file": self.preferred_cal_file.currentText()
        }
