"""
Design of the GUI for project creation
"""
from typing import Any, Dict

from PyQt5.QtWidgets import QLabel, QDialog, QLineEdit, QVBoxLayout, QDialogButtonBox, QComboBox


class ProjectSetupDialog(QDialog):
    """
    Design of the dialog box that appears for project creation. User can enter details in it.
    """

    def __init__(self, parent=None, interface=None):
        """
        Constructor

        :param parent: Parent
        :type parent: object
        :param interface: Tells which instrument specific interface is used. Eg: 'Metronix'
        :type interface: str
        """
        super().__init__(parent)
        self.interface = interface
        self.setWindowTitle("Project Setup")
        layout = QVBoxLayout()

        self.project_name = QLineEdit('XYZ Project')
        self.country = QLineEdit('India')
        self.acquired_by_company = QLineEdit('XYZ Organization')
        #
        self.processing_mode = QComboBox()
        self.processing_mode.addItem('MT + Tipper')
        self.processing_mode.addItem('MT Only')
        # TODO: Support use of tipper only in future
        # self.processing_mode.addItem('Tipper Only')
        #
        self.frequencies_per_decade = QLineEdit('12')
        self.notch_frequency = QLineEdit('50')
        self.preferred_cal_file = QComboBox()
        self.preferred_cal_file.addItem('xml')
        # TODO: Support use of metronix text file in future
        # self.preferred_cal_file.addItem('metronix txt')

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

    def get_values(self) -> Dict[str, Any]:
        """
        Return the values entered in the dialog box by the user.

        :return: Dictionary of content from the dialog box
        """
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
