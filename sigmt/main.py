"""
This is the first landing window for the SigMT
"""

import sys

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QPushButton, QDialog, QVBoxLayout, QLabel

from sigmt.gui.metronix.metronix import MainWindow as Metronix

if hasattr(Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)


class Welcome(QDialog):
    """
    Design of the welcome page. User can select the instrument based processing.
    """

    def __init__(self) -> None:
        """
        Constructor
        """
        super().__init__()
        self.metronix_button = None
        self.main_window = None
        self.init_ui()

    def init_ui(self) -> None:
        """
        Initial interface design

        :return: None
        """
        self.setWindowTitle('Welcome to SigMT | A Tool for Magnetotelluric Data Processing')
        self.setWindowIcon(QIcon(r'sigmt\images\sigmt.ico'))
        self.setGeometry(100, 100, 500, 500)
        layout = QVBoxLayout()

        layout.addStretch()
        # Add a welcome message label
        welcome_label = QLabel("Welcome to <b>SigMT</b>!", self)
        layout.addWidget(welcome_label, alignment=Qt.AlignHCenter)

        layout.addSpacing(50)

        select_label = QLabel("Select Instrument:", self)
        layout.addWidget(select_label, alignment=Qt.AlignHCenter)

        self.metronix_button = QPushButton('Metronix', self)
        layout.addWidget(self.metronix_button, alignment=Qt.AlignHCenter)

        layout.addStretch()

        self.setLayout(layout)

        self.metronix_button.clicked.connect(self.open_metronix_mainwindow)

    def open_metronix_mainwindow(self) -> None:
        """
        To open the metronix specific main window and close welcome page

        :return: None
        """
        self.main_window = Metronix()
        self.main_window.show()
        self.close()


def main() -> None:
    """
    Main Function

    :return: None
    """
    app = QApplication(sys.argv)
    _ = app.font()
    entry = Welcome()
    entry.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
