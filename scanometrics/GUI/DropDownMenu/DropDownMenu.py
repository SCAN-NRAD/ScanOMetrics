from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QComboBox

from scanometrics.GUI.Icon import Icon


class DropDownMenu(QWidget):
    """
    A custom drop-down menu widget with signal emission on option selection.
    """

    sgn = pyqtSignal(str)  # Signal emitted when an option is selected from the dropdown menu

    def __init__(self, default_msg, help_text):
        """
        Initialize the DropDownMenu widget.

        Args:
            default_msg (str): The default message to display in the dropdown menu.
            help_text (str): The tooltip text for the help icon.
        """
        super().__init__()
        self.layout = QHBoxLayout(self)  # Create the main horizontal layout

        self.help_icon = Icon(help_text)  # Create the help icon with customizable tooltip
        self.default_msg = default_msg  # Default message for the dropdown menu
        self.dropdown_menu = QComboBox(self)  # QComboBox instance for the dropdown menu
        self.initUI()  # Initialize the UI components

    def initUI(self):
        """
        Set up the dropdown menu within the layout and connect its activation signal.
        """
        self.set_default_message(self.default_msg)  # Set the initial default message

        # Add help icon and dropdown menu to the horizontal layout
        self.layout.addWidget(self.help_icon, stretch=2)
        self.layout.addWidget(self.dropdown_menu, stretch=98)

        self.setLayout(self.layout)  # Set the main layout for this widget

        self.dropdown_menu.activated.connect(self.select_option)  # Connect activation signal to handler

    def set_default_message(self, message):
        """
        Set and display the default message in the dropdown menu.

        Args:
            message (str): The message to be set as the default in the dropdown menu.
        """
        self.dropdown_menu.clear()  # Clear existing items
        self.dropdown_menu.addItem(message)  # Add the default message to the dropdown menu
        self.dropdown_menu.setCurrentText(message)  # Set the current text to the default message

    def select_option(self):
        """
        Handle the selection of an option from the dropdown menu.

        Emits the selected option as a signal if it is not the default message.
        """
        option = self.dropdown_menu.currentText()  # Get the currently selected option
        if option != self.default_msg:
            self.sgn.emit(option)  # Emit the selected option as a signal
            # print(option)  # Print the selected option for debugging
        else:
            self.sgn.emit(self.default_msg)  # Emit the default message as a signal
