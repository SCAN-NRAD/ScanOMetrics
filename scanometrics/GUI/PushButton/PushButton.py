from PyQt5.QtWidgets import QPushButton, QApplication, QStyle
from PyQt5.QtCore import pyqtSignal


class PushButton(QPushButton):
    """
    A custom push button widget with signal emission on click.
    """

    sgn = pyqtSignal(bool)  # Signal emitted when the button is clicked

    def __init__(self, label=""):
        """
        Initialize the PushButton widget.

        Args:
            label (str): The label text or key for predefined icons.
        """
        icon_map = {
            "left": QApplication.style().standardIcon(QStyle.SP_ArrowLeft),
            "right": QApplication.style().standardIcon(QStyle.SP_ArrowRight)
        }

        if label in icon_map:
            # If the label corresponds to a predefined icon, create the button with the icon
            super().__init__(icon_map[label], "")
        else:
            # Otherwise, create the button with the provided label text
            super().__init__(label)

        self.clicked.connect(self.on_click)  # Connect the button click event to the on_click method

    def on_click(self):
        """
        Emit the signal when the button is clicked.
        """
        self.sgn.emit(True)
