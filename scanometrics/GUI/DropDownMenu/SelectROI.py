from .DropDownMenu import DropDownMenu


class SelectROI(DropDownMenu):
    """
    Class for selecting Region of Interest (ROI) from a dropdown menu.
    """

    def __init__(self, default_msg):
        """
        Initialize the SelectROI class.

        Args:
            default_msg (str): The default message to display in the dropdown menu.
        """
        super().__init__(default_msg, help_text='This is the list with all the metrics depending on the checkbox clicked')

    def initUI(self):
        """
        Initialize the UI for the dropdown menu.
        """
        super().initUI()

    def get_items(self, metric_names):
        """
        Populate the dropdown menu with items based on the contents of a specified directory.

        Args:
            metric_names (list): A list of metric names to be added to the dropdown menu.
        """
        self.dropdown_menu.clear()
        self.dropdown_menu.addItem(self.default_msg)
        for name in metric_names:
            if 'aparc' not in name and 'symmetryIndex' not in name:
                self.dropdown_menu.addItem(name)

