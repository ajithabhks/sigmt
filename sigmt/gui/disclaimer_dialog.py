"""
Dialog box with disclaimer information (no HTML)
"""
import webbrowser

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QWidget, QHBoxLayout
)

from sigmt.__version__ import __version__


class DisclaimerDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.get_window_details()

        root = QVBoxLayout(self)

        # Version label
        version_label = QLabel(f"SigMT version: {__version__}")
        version_label.setAlignment(Qt.AlignCenter)
        version_label.setStyleSheet("font-weight:600;")

        # Warning "card"
        card = QWidget()
        card.setStyleSheet("""
            QWidget {
                background: #fff5f5;
                border: 1px solid #f5c2c7;
                border-left: 6px solid #dc3545;
                border-radius: 12px;
            }
        """)
        card_layout = QVBoxLayout(card)
        card_layout.setContentsMargins(16, 14, 16, 14)
        card_layout.setSpacing(10)

        # Title
        title = QLabel("⚠️ Important Disclaimer")
        title.setStyleSheet("font-size: 14px; font-weight: 700; color: #842029;")
        title.setAlignment(Qt.AlignLeft)

        card_layout.addWidget(title)

        # Bullet items (no HTML)
        bullets = [
            "This is the first experimental release of Phoenix Data Processing in SigMT.",
            "This feature has been tested only with MTU-5C data.",
            "The application may contain UI glitches and incorrect or incomplete results.",
            "Remote references have been tested, but may still produce unexpected issues.",
            "Report issues/suggestions via email (in GitHub page) or using the GitHub Issues page.",
            "This feature is under active development. A more stable release will be published after feedback.",
        ]

        for b in bullets:
            lbl = QLabel(f"•  {b}")
            lbl.setWordWrap(True)
            lbl.setStyleSheet("color: #842029;")
            card_layout.addWidget(lbl)

        # Links row (no HTML, clickable via mouse events)
        links_row = QHBoxLayout()
        links_row.setSpacing(12)

        github_link = self._make_link_label("SigMT GitHub Page", "https://github.com/ajithabhks/sigmt")
        issues_link = self._make_link_label("Report an Issue", "https://github.com/ajithabhks/sigmt/issues")

        links_row.addWidget(github_link)
        links_row.addWidget(issues_link)
        links_row.addStretch(1)

        card_layout.addLayout(links_row)

        # Layout
        root.addWidget(version_label)
        root.addWidget(card)

    def get_window_details(self) -> None:
        self.setWindowTitle("Important Disclaimer")
        self.setFixedSize(560, 420)
        self.setWindowIcon(QIcon(r"sigmt\images\sigmt.ico"))

    def _make_link_label(self, text: str, url: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setCursor(Qt.PointingHandCursor)
        lbl.setStyleSheet("""
            QLabel {
                color: #0b5ed7;
                text-decoration: underline;
                font-weight: 600;
            }
            QLabel:hover {
                color: #084298;
            }
        """)

        # Attach a click handler without HTML
        def open_link(_event):
            webbrowser.open(url)

        lbl.mousePressEvent = open_link
        return lbl
