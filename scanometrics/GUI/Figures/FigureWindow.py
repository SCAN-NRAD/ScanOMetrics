import io
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QCheckBox, QMessageBox
from PyQt5.QtCore import pyqtSignal, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from .SpiderPlot import *


class FigureWindow(QDialog):
    sgn = pyqtSignal(object)  # Signal to emit the buffed figure

    def __init__(self, parent=None, fig_type=None):
        super().__init__(parent)
        self.parent, self.fig_type = parent, fig_type
        self.main_layout, self.checkbox_layout = QVBoxLayout(), QHBoxLayout()
        self.fig_cloud = self.fig_spider = self.canvas_cloud = self.canvas_spider = None
        self.checkbox_cloud, self.checkbox_spider = QCheckBox('Add Cloud Point'), QCheckBox('Add Spider Plot')

        # Configure and set the figures
        self.set_figures()

        # Build the figures window
        self.window_builder()

        self.setLayout(self.main_layout)

        # Show the dialog window
        self.show()

    def set_figures(self):
        """Sets the cloud point and spider plot figures."""
        self.setWindowTitle(f"{self.parent.sub_name} - Figures")
        if self.fig_type == 'cloud_points':
            self.get_cloud_points()
        else:
            if self.parent.check_value:
                aparc_code = 'aparc' if self.parent.brain3d.atlas == 'DesikanKilliany' else 'aparca2009s'
                self.parent.metric = f"{aparc_code}_lh_{self.parent.brain_region_info['region_name']}_{self.parent.check_value}"
                self.get_cloud_points()
                self.get_spiderplot()
            else:
                self.get_spiderplot()

    def get_cloud_points(self):
        """
        Create a cloud point figure
        """
        self.fig_cloud = self.parent.som_project.plot_single_metric(
            self.parent.sub_name, self.parent.metric, self.parent.ses_nb, self.parent.acq_lab,
            self.parent.evaluation_results['output'],
            # self.parent.evaluation_results['subj_covariate_values'],
            # self.parent.evaluation_results['normModel_matching_cols'],
            # self.parent.evaluation_results['subj_measured_metrics'],
            ['sex', 'sequence'],  # TODO: add option to set covariates in GUI (and sync with call to evaluation)
            self.parent.evaluation_results['pid'], 0.01)

    def get_spiderplot(self):
        """
        Create a spider plot figure
        """
        metrics, lh_stats, rh_stats, asymmetry_stats = get_values(self.parent.brain_region_info,
                                                                  self.parent.evaluation_results['output'])
        self.fig_spider, axs = plt.subplots(1, 3, figsize=(12, 6), subplot_kw=dict(polar=True))
        self.fig_spider.suptitle(f"{self.parent.brain_region_info['region_name']}", color='black', fontsize=16)
        create_spiderplot(axs[0], metrics, lh_stats, "Left Hemisphere")
        create_spiderplot(axs[1], metrics, asymmetry_stats, "Asymmetry")
        create_spiderplot(axs[2], metrics, rh_stats, "Right Hemisphere")
        handles, labels = axs[0].get_legend_handles_labels()
        self.fig_spider.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.95, 0.05))

    def window_builder(self):
        # Set up the dialog window
        self.setWindowTitle('Figures')
        self.setGeometry(200, 200, 1000, 800)
        self.setWindowFlags(
            Qt.Window | Qt.WindowMaximizeButtonHint | Qt.WindowMinimizeButtonHint | Qt.WindowCloseButtonHint)

        [self.main_layout.addWidget(FigureCanvas(fig)) for fig in [self.fig_cloud, self.fig_spider] if fig]
        [self.checkbox_layout.addWidget(checkbox) for fig, checkbox
         in zip([self.fig_cloud, self.fig_spider], [self.checkbox_cloud, self.checkbox_spider]) if fig]
        self.main_layout.addLayout(self.checkbox_layout)
        self.checkbox_cloud.setChecked(True)
        self.checkbox_spider.setChecked(True)
        self.main_layout.addWidget(button := QPushButton('Add to the report'))
        button.clicked.connect(self.emit_fig_data)

    def emit_fig_data(self):
        """Emits the selected figure data as a signal and closes the dialog."""
        selected_figs = []

        if self.checkbox_cloud.isChecked() and self.fig_cloud:
            buf_cloud = io.BytesIO()
            self.fig_cloud.savefig(buf_cloud, format='png')
            buf_cloud.seek(0)
            selected_figs.append({'image': buf_cloud, 'caption': f"{self.parent.sub_name}_cloud_point"})

        if self.checkbox_spider.isChecked() and self.fig_spider:
            buf_spider = io.BytesIO()
            self.fig_spider.savefig(buf_spider, format='png')
            buf_spider.seek(0)
            selected_figs.append({'image': buf_spider, 'caption': f"{self.parent.sub_name}_spiderplot"})

        if selected_figs:
            self.sgn.emit(selected_figs)
            self.clear_layout(self.main_layout)
            self.canvas_cloud = self.canvas_spider = None
            self.fig_cloud = self.fig_spider = None
            self.close()
        else:
            QMessageBox.warning(self, 'No Selection', 'Please select at least one figure to add to the report.')

    def clear_layout(self, layout):
        """Clear all widgets from the given layout."""
        if layout is not None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clear_layout(child.layout())
