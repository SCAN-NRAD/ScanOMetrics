import os
from .DropDownMenu import DropDownMenu


class SelectSub(DropDownMenu):
    """
    A dropdown menu widget specifically designed for selecting subject directories.
    """

    def __init__(self, default_msg):
        """
        Initialize the SelectSub widget.

        Args:
            default_msg (str): The default message to display in the dropdown menu.
        """
        super().__init__(default_msg, help_text='Here you need to choose one subject for the evaluation')

    def initUI(self):
        """
        Set up the UI by calling the parent initUI method.
        """
        super().initUI()

    def get_items(self, sample_folder_path):
        """
        Populate the dropdown menu with subject directories.

        Args:
            sample_folder_path (str): The path to the sample folder.
        """
        self.dropdown_menu.clear()
        self.dropdown_menu.addItem(self.default_msg)
        try:
            subjects = [sub for sub in os.listdir(sample_folder_path) if os.path.isdir(os.path.join(sample_folder_path, sub)) and "sub" in sub]
            subjects.sort()
            if len(subjects) == 1:
                self.dropdown_menu.addItem(subjects[0])
                self.dropdown_menu.setCurrentIndex(1)
                self.dropdown_menu.activated.emit(1)
            else:
                for sub in subjects:
                    self.dropdown_menu.addItem(sub)
        except Exception as e:
            print(f"Error accessing sample folder: {sample_folder_path} {e}")
