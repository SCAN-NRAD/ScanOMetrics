from PyQt5.QtWidgets import QWidget, QVBoxLayout, QCheckBox
from PyQt5.QtCore import pyqtSignal, Qt


class CheckBox(QWidget):
    # Signal emitted when the checkbox state changes. It sends the label of the checkbox when checked.
    sgn = pyqtSignal(str)

    def __init__(self, label):
        super().__init__()
        # Initialize the QCheckBox with the provided label
        self.checkbox = QCheckBox(label)

        # Connect the checkbox's stateChanged signal to the custom slot emit_value_changed
        self.checkbox.stateChanged.connect(self.emit_signal)

        # Set up the layout for this custom widget
        layout = QVBoxLayout()
        layout.addWidget(self.checkbox)
        self.setLayout(layout)

    def emit_signal(self, state):
        """
        Slot that gets called when the checkbox state changes.
        Emits the valueChanged signal with the checkbox's label only if the checkbox is checked.
        """
        if state == Qt.Checked:
            # Emit the signal with the checkbox's text when it is checked
            self.sgn.emit(self.checkbox.text())

    def set_checked(self, checked):
        """
        Method to programmatically set the checkbox state.
        :param checked: Boolean indicating whether the checkbox should be checked (True) or not (False).
        """
        self.checkbox.setChecked(checked)

    def is_checked(self):
        """
        Method to check the current state of the checkbox.
        :return: True if the checkbox is checked, otherwise False.
        """
        return self.checkbox.isChecked()
