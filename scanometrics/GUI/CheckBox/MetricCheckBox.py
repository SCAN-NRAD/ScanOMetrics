from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget, QHBoxLayout

from .CheckBox import CheckBox


class MetricCheckBox(QWidget):
    # Signal emitted when the checkbox state changes. It sends the label of the checkbox when checked.
    sgn = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.checkboxes_layout = QHBoxLayout()
        self.selected_value = None  # To store the value of the selected checkbox

        widgets_info = [
            ('volume_box', CheckBox('volume'), 'color: #ffffff', 'volume', 'left'),
            ('thickness_box', CheckBox('thickness'), 'color: #ffffff', 'thickness', 'left'),
            ('thicknessstd_box', CheckBox('thicknessstd'), 'color: #ffffff', 'thicknessstd', 'left')
        ]
        self.fill_layout(widgets_info)

        self.setLayout(self.checkboxes_layout)

    def fill_layout(self, widgets_info):
        """Helper method to fill the specified layout with widgets.

        Args:
            widgets_info (list): List of tuples containing widget information:
                                 (widget_name, widget, widget_style, signal_label, layout_position)
        """
        # Iterate through widget information and add widgets to the corresponding layout
        for widget_name, widget, widget_style, sgn_label, pos in widgets_info:
            widget.setStyleSheet(widget_style)
            if sgn_label:
                # Connect the widget's signal to the handle_sgn slot
                widget.sgn.connect(self.handle_sgn)
            self.checkboxes_layout.addWidget(widget)
            setattr(self, widget_name, widget)  # Store widget as an attribute of MainWindow

    def handle_sgn(self, sgn_value):
        """Handle the signal emitted by the checkboxes and store the selected value.

        Args:
            sgn_value (str): The value of the checkbox that was checked.
        """
        self.selected_value = sgn_value
        # Emit the signal from MetricCheckBox
        self.sgn.emit(self.selected_value)

    def reset_checkboxes(self):
        """Reset all checkboxes except the one with the current value."""
        for widget_name in ['volume_box', 'thickness_box', 'thicknessstd_box']:
            checkbox = getattr(self, widget_name)
            if checkbox.checkbox.text() != self.selected_value:
                checkbox.set_checked(False)

    def get_selected_value(self):
        """Return the currently selected checkbox value."""
        return self.selected_value
