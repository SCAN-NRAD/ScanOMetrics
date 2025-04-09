from PyQt5.QtWidgets import QLabel, QApplication, QStyle, QToolTip
from PyQt5.QtCore import Qt


class Icon(QLabel):
    """
    A custom widget representing an icon with a tooltip.
    """

    def __init__(self, tooltip_text, icon_name='SP_MessageBoxQuestion'):
        """
        Initialize the Icon widget.

        Args:
            tooltip_text (str): The tooltip text for the icon.
            icon_name (str): The name of the standard icon to use.
        """
        super().__init__()
        self.tooltip_text = tooltip_text  # Store the tooltip text

        try:
            style = QApplication.style()  # Get the application style
            icon = style.standardIcon(getattr(QStyle, icon_name))  # Get a standard icon
            if icon.isNull():
                raise FileNotFoundError("Standard icon not found.")  # Raise an exception if the icon is not found
            self.setPixmap(icon.pixmap(20, 20))  # Set the icon pixmap with a specific size
            self.setFixedSize(20, 20)  # Set the QLabel size to match the icon size
        except Exception as e:
            print(f"Error loading icon: {e}")  # Print an error message if there's an issue loading the icon

    def enterEvent(self, event):
        """
        Event handler for mouse entering the widget.

        Args:
            event: The event object.
        """
        QToolTip.showText(event.globalPos(), self.tooltip_text, self)
        super().enterEvent(event)