import importlib
import os
from PyQt5.QtWidgets import QMessageBox
from .DropDownMenu import DropDownMenu


class SelectPKLFile(DropDownMenu):
    """
    A dropdown menu widget specifically designed for selecting .pkl files.
    """

    def __init__(self, default_msg, normative_model_ids):
        """
        Initialize the SelectPKLFile widget.

        Args:
            default_msg (str): The default message to display in the dropdown menu.
        """
        super().__init__(default_msg, help_text="The .pkl file is the normative model")
        self.add_items(normative_model_ids)

    def initUI(self):
        """
        Set up the UI and launch the .pkl file search.
        """
        super().initUI()  # Call the parent initUI method to set up the basic UI

    def add_items(self, normative_model_ids):
        """
        Locate .pkl files in the specified directory and add them to the dropdown menu.
        """
        self.dropdown_menu.clear()
        self.dropdown_menu.addItem(self.default_msg)
        try:
            if len(normative_model_ids) == 1:
                self.dropdown_menu.addItem(normative_model_ids[0])
                self.dropdown_menu.setCurrentText(normative_model_ids[0])
                self.sgn.emit(normative_model_ids[0])
            else:
                for normative_model_id in normative_model_ids:
                    self.dropdown_menu.addItem(normative_model_id)
        except Exception as e:
            # If there's an error accessing the directory, display an error message
            QMessageBox.critical(self, "Error", f"Error while setting list of normative models: {e}")
