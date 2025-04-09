import os
import re
import traceback

from glob import glob
import io
import vtk
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import pyqtSignal, QBuffer
from PyQt5.QtWidgets import QHBoxLayout, QLabel, QVBoxLayout, QPushButton, QMessageBox
import nibabel.freesurfer.io as fs
from scanometrics.utils import logging
from scanometrics.GUI import styles
from scanometrics.GUI.Brain3D.Hemisphere import Hemisphere
from scanometrics.GUI.PushButton import PushButton


class ColorBar(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setFixedWidth(25)
        self.setFixedHeight(375)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        gradient = QtGui.QLinearGradient(0, 0, 0, self.height())

        gradient.setColorAt(0.0, QtGui.QColor(255, 255, 0))  # Yellow
        gradient.setColorAt(0.2, QtGui.QColor(255, 0, 0))  # Red
        gradient.setColorAt(0.215, QtGui.QColor(80, 80, 80))  # Dark Gray
        gradient.setColorAt(0.8225, QtGui.QColor(80, 80, 80))  # Dark Gray
        gradient.setColorAt(0.825, QtGui.QColor(0, 0, 255))  # Blue
        gradient.setColorAt(1.0, QtGui.QColor(0, 255, 255))  # Light Blue

        painter.fillRect(self.rect(), gradient)

        painter.setPen(QtGui.QPen(QtCore.Qt.black, 2))
        metrics = painter.fontMetrics()
        labels = [2.0, 1.3010299956639813, 0, -1.3010299956639813, -2.0]
        positions = [0.0, 0.2, 0.5, 0.8, 1.0]

        for i, label in enumerate(labels):
            text = f"{label:.1f}"
            text_width = metrics.width(text)
            painter.drawText(int((self.width() - text_width) / 2), int(positions[i] * self.height() + metrics.ascent()),
                             text)

aparc_code = {'DesikanKilliany': 'aparc',
              'Destrieux': 'aparca2009s'}

class Brain3D(QtWidgets.QWidget):
    sgn = pyqtSignal(dict)

    def __init__(self, scan_dir):
        """
        Create Brain3D rendering for a given subject.
        The subject directory (TODO: define directory variable from MainWindow) and evaluation output (self.output) are set from MainWindow
        Evaluation results are passed from the EvaluationThread into
        """
        super().__init__()
        self.output = self.lh_vals = self.rh_vals = None
        self.brain_model = self.brain_annot = ()
        self.brain_layout, self.button_layout, self.hemispheres_layout = QVBoxLayout(), QHBoxLayout(), QHBoxLayout()
        self.sync_enabled = False  # Attribute to track the hemisphere's camera synchronization state
        # DL+DiReCT saves both DesikanKilliany and Destrieux labels as aparc.annot, while Freesurfer saves them as aparc.annot and aparc.a2009s.annot respectively
        if 'dldirect' in scan_dir:
            self.pipeline = 'dldirect'
        elif 'freesurfer' in scan_dir:
            self.pipeline = 'freesurfer'
        if 'DesikanKilliany' in scan_dir:
            self.atlas = 'DesikanKilliany'
        elif 'Destrieux' in scan_dir:
            self.atlas = 'Destrieux'
        if os.path.exists(os.path.join(scan_dir, 'surf', 'lh.pial')) and os.path.exists(os.path.join(scan_dir, 'surf', 'rh.pial')) and os.path.exists(os.path.join(scan_dir, 'label', 'lh.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s'))) and os.path.exists(os.path.join(scan_dir, 'label', 'rh.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s'))):
            self.subject_dir = scan_dir
            self.brain_model = tuple(os.path.join(self.subject_dir, 'surf', f'{hemi}.pial') for hemi in ['lh', 'rh'])
            self.brain_annot = tuple(os.path.join(self.subject_dir, 'label', f'{hemi}.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s')) for hemi in ['lh', 'rh'])
        else:
            if not all(os.path.exists(os.path.join(scan_dir, 'surf', hemi + '.pial')) for hemi in ['lh', 'rh']):
                logging.WARNING('No pial surface found in %s/surf, ' % scan_dir)
            if not os.path.exists(os.path.join(scan_dir, 'label', 'lh.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s'))) or not os.path.exists(os.path.join(scan_dir, 'label', 'rh.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s'))):
                logging.WARNING('No annotation file found in %s/label' % scan_dir)
            if "FREESURFER_HOME" in os.environ:
                self.subject_dir = os.path.join(os.environ["FREESURFER_HOME"], 'subjects', 'fsaverage')
                logging.WARNING('Falling back to FREESURFER fsaverage rendering (%s)' % self.subject_dir)
                self.brain_model = tuple(os.path.join(self.subject_dir, 'surf', f'{hemi}.pial') for hemi in ['lh', 'rh'])
                self.brain_annot = tuple(os.path.join(self.subject_dir, 'label', f'{hemi}.%s.annot' % ('aparc' if self.atlas == 'DesikanKilliany' else 'aparc.a2009s')) for hemi in ['lh', 'rh'])
            else:
                logging.ERROR("%s is missing surface/label files and FREESURFER_HOME not found. Please set it up before running ScanOMetrics (eg: 'export FREESURFER_HOME=/usr/local/freesurfer;source /usr/local/freesurfer/SetUpFreeSurfer.sh')")

        widgets_info = [
            ('sync_button', PushButton("Sync Cameras"), styles.button_style, 'synchro', 'button'),
            ('reset_button', PushButton("Reset Cameras"), styles.button_style, 'reset', 'button'),
            ('left_brain', Hemisphere(self.brain_model[0],
                                      self.brain_annot[0],
                                      [(-512, 128, 128),
                                                 (0, 0, 0),
                                                 (0, 0, 1)]),
                           None, 'click_brain', 'hemisphere'),
            ('color_bar', ColorBar(), None, None, 'hemisphere'),
            ('right_brain', Hemisphere(self.brain_model[1],
                                       self.brain_annot[1],
                                       [(512, 128, 128),
                                                  (0, 0, 0),
                                                  (0, 0, 1)]),
                            None, 'click_brain', 'hemisphere'),
            ('region_label', QLabel("Click on a brain region to view information"), styles.brain_label_style, None, 'brain'),
            ('save_button', PushButton("Capture the brain"), styles.button_style, 'save', 'brain')
        ]
        self.fill_layout(widgets_info)
        self.setLayout(self.brain_layout)
        self.sync_cam()

    def fill_layout(self, widgets_info):
        """Method to fill the specified layout with widgets.

        Args:
            widgets_info (list): List of tuples containing widget information:
                                 (widget_name, widget, widget_style, signal_label, layout_position)
        """
        # Mapping layout positions to actual layout objects
        pos_lay = {
            'brain': self.brain_layout,
            'button': self.button_layout,
            'hemisphere': self.hemispheres_layout
        }

        # Iterate through widget information and add widgets to the corresponding layout
        for widget_name, widget, widget_style, sgn_label, pos in widgets_info:
            if widget_style:
                widget.setStyleSheet(widget_style)
                if widget_style == styles.brain_label_style:
                    widget.setAlignment(QtCore.Qt.AlignCenter)
            if sgn_label:
                # Connect the widget's signal to the handle_sgn slot, passing the signal label
                widget.sgn.connect(lambda sgn, lbl=sgn_label: self.handle_sgn(sgn, lbl))
            pos_lay[pos].addWidget(widget)
            setattr(self, widget_name, widget)  # Store widget as an attribute of MainWindow
            if widget_name == 'reset_button':
                self.brain_layout.addLayout(self.button_layout)
            if widget_name == 'right_brain':
                self.brain_layout.addLayout(self.hemispheres_layout)

    def handle_sgn(self, value, signal_type):
        try:
            if signal_type == 'synchro':
                self.sync_cam()
            if signal_type == 'reset':
                self.reset_cam()
            if signal_type == 'click_brain':
                self.update_region_label(value)
            if signal_type == 'save':
                self.capture_and_emit_screenshot()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}\n\n{traceback.format_exc()}")

    def sync_cam(self):
        """Synchronize or unsynchronize the camera movements between the two hemispheres."""
        self.sync_enabled = not self.sync_enabled
        if self.sync_enabled:
            self.sync_button.setStyleSheet(styles.button_checked)
            self.left_brain.camera_updated.connect(self.right_brain.update_camera)
            self.right_brain.camera_updated.connect(self.left_brain.update_camera)
        else:
            self.sync_button.setStyleSheet(styles.button_style)
            self.left_brain.camera_updated.disconnect(self.right_brain.update_camera)
            self.right_brain.camera_updated.disconnect(self.left_brain.update_camera)

    def reset_cam(self):
        """Reset cameras to their initial positions."""
        self.left_brain.reset_camera_position()
        self.right_brain.reset_camera_position()

    def update_region_label(self, region_info):
        self.region_label.setText(f"Selected Region: {region_info['region_name']},"
                                  f" Coordinates: {region_info['picked_position']}")
        self.sgn.emit(region_info)  # Emit signal with region info

    def update_display(self, metric):
        """
        Update the 3D brain model display with the latest annotation colors.
        """
        self.get_logps_vals(metric)
        for hemisphere, dic_devtn_logps in {self.left_brain: self.lh_vals, self.right_brain: self.rh_vals}.items():
            hemisphere.update_polydata_with_annot(dic_devtn_logps)
            hemisphere.vtkWidget.GetRenderWindow().Render()

    def get_logps_vals(self, metric):
        labels, ctab, region_names = fs.read_annot(self.brain_annot[0])
        region_names = [region.decode('utf-8') for region in region_names]
        for prefix in ['lh', 'rh']:
            devtn_logps = []
            for region in region_names:
                metric_name = f"%s_{prefix}_{region}_{metric}" % (aparc_code[self.atlas])
                if metric_name in self.output['norm']['metric_names']:
                    logp_value = self.output['norm']['devtn_logps'][0, self.output['norm']['metric_names'].index(metric_name)]
                    devtn_logps.append(logp_value)
                else:
                    devtn_logps.append(0)
            setattr(self, f"{prefix}_vals", devtn_logps)

    def capture_and_emit_screenshot(self):
        """
        Capture screenshots from both hemispheres and emit them.
        """
        left_screenshot = self.left_brain.save_screenshot()
        right_screenshot = self.right_brain.save_screenshot()
        pixmap = self.color_bar.grab()
        buffer = QBuffer()
        buffer.open(QBuffer.ReadWrite)
        pixmap.save(buffer, "PNG")
        colorbar_image = io.BytesIO(buffer.data())
        buffer.close()

        screenshots = {
            'type': 'screenshots',
            'left_hemisphere': left_screenshot,
            'right_hemisphere': right_screenshot,
            'color_bar': {'image': colorbar_image}
        }

        self.sgn.emit(screenshots)
