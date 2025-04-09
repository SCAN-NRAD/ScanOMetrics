import os
from .DropDownMenu import DropDownMenu


class SelectSes(DropDownMenu):
    """
    A dropdown menu widget specifically designed for selecting session directories.
    """

    def __init__(self, default_msg):
        """
        Initialize the SelectSes widget.

        Args:
            default_msg (str): The default message to display in the dropdown menu.
        """
        super().__init__(default_msg, help_text="Here you need to choose the participant's session")

    def initUI(self):
        """
        Set up the UI by calling the parent initUI method.
        """
        super().initUI()

    def get_items(self, sample_folder_path, sub_name):
        """
        Populate the dropdown menu with session directories.

        Args:
            sample_folder_path (str): The path to the sample folder.
            sub_name (str): The name of the subject.
        """
        self.dropdown_menu.clear()
        self.dropdown_menu.addItem(self.default_msg)
        try:
            sessions_path = os.path.join(sample_folder_path, sub_name)
            sessions = [ses for ses in os.listdir(sessions_path) if os.path.isdir(os.path.join(sessions_path, ses)) and "ses" in ses]
            sessions.sort()
            if len(sessions) == 1:
                self.dropdown_menu.addItem(sessions[0])
                self.dropdown_menu.setCurrentIndex(1)
                self.dropdown_menu.activated.emit(1)
            else:
                for ses in sessions:
                    self.dropdown_menu.addItem(ses)
        except Exception as e:
            print(f"Error accessing {sub_name} folder - {e}")
