import os
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QPushButton, QFileDialog, QMessageBox
from PyQt5.QtCore import pyqtSignal
from scanometrics.GUI import styles
from scanometrics.GUI.Icon import Icon


class PathSelector(QWidget):
    """
    A custom widget for selecting a folder path with additional validation.
    """

    sgn = pyqtSignal(str)  # Signal emitted when a path is selected

    def __init__(self, placeholder_text, help_text="The sample folder's format must be in OASIS like format"):
        """
        Initialize the PathSelector widget.

        Args:
            placeholder_text (str): The placeholder text for the path line edit.
            help_text (str): The tooltip text for the help icon.
        """
        super().__init__()
        self.layout = QHBoxLayout(self)  # Create the main horizontal layout
        self.help_icon = Icon(help_text)  # Create the help icon with customizable tooltip
        self.path_line = QLineEdit(self)  # Create the line to write the folder's path
        self.browse_button = QPushButton("Search", self)  # Create the browse button to search the folder's path
        self.initUI(placeholder_text)

    def initUI(self, placeholder_text):
        """
        Initialize the user interface.

        Args:
            placeholder_text (str): The placeholder text for the path line edit.
        """
        self.path_line.setPlaceholderText(placeholder_text)  # Set placeholder text for the path line edit

        self.browse_button.setStyleSheet(styles.button_style)  # Apply the button style

        self.browse_button.clicked.connect(self.browse_path)  # Connect browse button click event to browse_path method

        # Add help icon, QLineEdit, and QPushButton to the horizontal layout
        self.layout.addWidget(self.help_icon, stretch=2)
        self.layout.addWidget(self.path_line, stretch=90)
        self.layout.addWidget(self.browse_button, stretch=8)

        self.setLayout(self.layout)  # Set the main layout for this widget

    def browse_path(self):
        """Open a directory selection dialog and emit the selected path."""
        try:
            # Open directory selection dialog and get the selected path
            path = QFileDialog.getExistingDirectory(self, "Select Sample Folder")
            if path:
                # Check the folder structure and emit the selected path if valid
                self.check_bids_format(path)
                self.path_line.setText(path)  # Display the selected path in the line edit
                self.sgn.emit(path)  # Emit the selected path signal
            else:
                # Display a warning message if no directory is selected
                QMessageBox.warning(self, "No Directory Selected", "Please select a valid directory.")
        except Exception as e:
            # Display an error message if an error occurs during directory selection
            QMessageBox.critical(self, "Error", f"An error occurred while selecting the directory:\n{e}")

    @staticmethod
    def check_bids_format(path):
        """
        Validate the structure of a directory to ensure it matches a bids-like structure.

        Args:
            path (str): The path to the directory to be checked.
        """
        # Check for the presence of required directories and files
        if not os.path.isdir(os.path.join(path, 'code')):
            raise Exception("The 'code' directory is missing.")
        if not os.path.isdir(os.path.join(path, 'derivatives')):
            raise Exception("The 'derivatives' directory is missing.")
        if not os.path.isfile(os.path.join(path, 'participants.tsv')):
            raise Exception("The 'participants.tsv' file is missing.")

        # Check for 'sub' directories and 'ses' directories within them
        sub_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and 'sub' in d]
        if not sub_dirs:
            raise Exception("No directory with 'sub' in the name was found.")
        for sub_dir in sub_dirs:
            ses_dirs = [d for d in os.listdir(os.path.join(path, sub_dir)) if
                        os.path.isdir(os.path.join(path, sub_dir, d)) and 'ses' in d]
            if not ses_dirs:
                raise Exception(f"No directory with 'ses' in the name was found in {sub_dir}.")
            for ses_dir in ses_dirs:
                anat_dir_path = os.path.join(path, sub_dir, ses_dir, 'anat')
                if not os.path.isdir(anat_dir_path):
                    raise Exception(f"The 'anat' directory is missing in {os.path.join(path, sub_dir, ses_dir)}.")
                if not os.listdir(anat_dir_path):
                    raise Exception(f"The 'anat' directory is empty in {os.path.join(path, sub_dir, ses_dir)}.")
