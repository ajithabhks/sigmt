"""
Dialog box with information about the package
"""
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel

from sigmt.__version__ import __version__


class AboutDialog(QDialog):
    """
    Dialog box with information about the package
    """

    def __init__(self):
        super().__init__()

        # Set the title and size of the dialog
        self.citation_label = None

        self.get_window_details()

        # Create a layout for the dialog
        layout = QVBoxLayout()

        self.version_label = QLabel(f'SigMT version: {__version__}')
        self.version_label.setAlignment(Qt.AlignCenter)

        self.get_citation_label()
        self.citation_label.setOpenExternalLinks(True)  # Enable clickable link
        self.citation_label.setAlignment(Qt.AlignCenter)

        # Add the label to the layout
        layout.addWidget(self.version_label)
        layout.addWidget(self.citation_label)

        # Set the layout for the dialog
        self.setLayout(layout)

    def get_window_details(self) -> None:
        """
        Getting window details
        """
        self.setWindowTitle("About")
        self.setFixedSize(400, 200)
        self.setWindowIcon(QIcon(r'sigmt\images\sigmt.ico'))

    def get_citation_label(self) -> None:
        """
        Create citation label
        """
        # Add the citation label with a clickable DOI link using rich text (HTML)
        self.citation_label = QLabel('''
                                     If you are using SigMT, please cite the following paper:<br><br>
                                     Ajithabh, K. S., Patro, P. K. (2023). SigMT: An open-source 
                                     Python package for magnetotelluric data processing. <br>
                                     <i>Computers & Geosciences</i>, 171, 105270.<br>
                                     DOI: <a href="https://doi.org/10.1016/j.cageo.2022.105270">
                                     https://doi.org/10.1016/j.cageo.2022.105270</a>
                                     ''')
