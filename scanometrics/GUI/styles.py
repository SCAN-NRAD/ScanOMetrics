background_style = """
QMainWindow {
    background: qlineargradient(x1: 0, y1: 0, x2: 1, y2: 1,
        stop: 0 #000000,
        stop: 0.25 #0f3460,
        stop: 0.5 #000000,
        stop: 0.75 #0f3460,
        stop: 1 #000000);
}
"""

separator_style = """
color: #4b77be; 
background-color: #4b77be;
"""

title_style = """
QLabel {
    color: #ffffff;
    font-size: 40px;
    font-weight: bold;
    font-family: Arial, sans-serif;
    text-align: center;
    padding: 20px;
    text-transform: uppercase;
    letter-spacing: 3px;
}
"""


button_style = """
QPushButton {
    background: transparent;
    color: #ffffff;
    border: 2px solid #4b77be;
    padding: 10px 20px;
    border-radius: 15px;
    font-size: 18px;
}
QPushButton:hover {background: rgba(255, 255, 255, 0.1);}
QPushButton:pressed {background: rgba(255, 255, 255, 0.2);}
"""
button_checked = """
QPushButton {
    background: rgba(255, 255, 255, 0.2);
    color: #ffffff;
    border: 2px solid #4b77be;
    padding: 10px 20px;
    border-radius: 15px;
    font-size: 18px;
}
QPushButton:hover {background: rgba(255, 255, 255, 0.1);}
QPushButton:pressed {background: rgba(255, 255, 255, 0.2);}
"""

line_edit_style = """
QLineEdit {
    background: #2c3e50;
    color: #ffffff;
    border: none;
    padding: 8px;
    border-radius: 15px;
    font-size: 16px;
}
"""

combobox_style = """
QComboBox {
    background: #2c3e50;
    color: #ffffff;
    border: none;
    padding: 8px;
    border-radius: 15px;
    font-size: 16px;
}
QComboBox QAbstractItemView {
    background-color: #1e1e1e;  
    color: #ffffff;
    selection-background-color: #4b77be;  
    selection-color: #ffffff;
    border: none;
    border-radius: 15px;
}
QComboBox::drop-down {
    border: none;
    background: transparent;
}
QComboBox QAbstractItemView::item {
    background: transparent;
    color: #ffffff;
}
QComboBox QAbstractItemView::item:hover {
    background: #34495e;
    color: #ffffff;
}
"""

brain_label_style = """
QLabel {
    color: #ffffff;
    font-size: 20px;
    text-align: center;
}
"""
