import os
import re
from .DropDownMenu import DropDownMenu


class SelectAcqLab(DropDownMenu):
    """
    A dropdown menu widget for selecting acquisition labels.
    """

    def __init__(self, default_msg):
        super().__init__(default_msg, help_text="Here you need to choose the acquisition label")

    def initUI(self):
        super().initUI()

    def get_items(self, sample_folder_path, sub_name, ses_nb):
        """
        Populate the dropdown menu with unique acquisition labels (T{x}w) from MRI files in the session directory.

        Args:
            sample_folder_path (str): The path to the sample folder.
            sub_name (str): The name of the subject.
            ses_nb (str): The session number.
        """
        self.dropdown_menu.clear()
        self.dropdown_menu.addItem(self.default_msg)
        ses_path = os.path.join(sample_folder_path, sub_name, ses_nb, "anat")

        if not os.path.exists(ses_path):
            print(f"Session directory does not exist: {ses_path}")
            return

        try:
            acq_parts = [re.search(r'^%s_%s_(.*).nii.gz' % (sub_name, ses_nb), file).group(1) for file in os.listdir(ses_path) if file.endswith(".nii.gz")]
            unique_acq_parts = list(set(acq_parts))
            unique_acq_parts.sort()
            if len(unique_acq_parts) == 1:
                self.dropdown_menu.addItem(unique_acq_parts[0])
                self.dropdown_menu.setCurrentIndex(1)
                self.dropdown_menu.activated.emit(1)
            else:
                for acq_part in unique_acq_parts:
                    self.dropdown_menu.addItem(acq_part)

        except Exception as e:
            print(f"An error occurred while processing files: {e}")
