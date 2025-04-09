import io
import sys
import traceback
from glob import glob
import os

from scanometrics.GUI.Brain3D.Brain3D import aparc_code
from scanometrics.utils.logging import set_verbose
from scanometrics.utils.zenodo_api import list_normative_models
import pickle


set_verbose(2)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel, QMessageBox, QGridLayout, \
    QHBoxLayout, QFrame, QProgressDialog
from PyQt5.QtCore import Qt, pyqtSlot

from scanometrics.GUI.Thread import EvaluationThread
from scanometrics.core import ScanOMetrics_project
from scanometrics.GUI.Brain3D import Brain3D
from scanometrics.GUI.PathSelector import PathSelector
from scanometrics.GUI.DropDownMenu import SelectPKLFile, SelectSub, SelectSes, SelectAcqLab, SelectROI
from scanometrics.GUI.PushButton import PushButton
from scanometrics.GUI.CheckBox import MetricCheckBox
from scanometrics.GUI.Figures import FigureWindow
from scanometrics.GUI.ReportGenerator import ReportGenerator
from scanometrics.GUI import styles


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # Initialize main UI components
        self.figure_dialog = None
        self.brain_region_info = None
        self.report_generator = None
        self.main_widget = self.main_layout = self.title_layout = self.title_label = \
            self.header_layout = self.left_layout = self.left_label = self.right_layout = \
            self.current_page = None

        # Initialize user input variables with default values
        self.normative_models = list_normative_models()
        (self.sample_folder_path, self.normative_model_string, self.sub_name,
         self.ses_nb, self.acq_lab, self.check_value, self.metric) = \
            ("Select the sample folder path...", "Select your normative model (.pkl)", "Select the Subject",
             "Select the Session", "Select acquisition label", "None", "Select the metrics")

        # Initialize ScanOMetrics project and 3D brain visualization instances
        self.som_project = self.brain3d = None

        # Initialize Thread components and the evaluate output variable
        self.evaluation_thread = self.progress_dialog = self.output = None

        # Initialize dictionary to store evaluation results
        self.evaluation_results = {}

        # List to store plot data
        self.plot_data_list = []

        # List to store brain screenshots
        self.brain_screenshots = []

        # Flag to check if report is created
        self.is_report_created = False

        # Set up the user interface
        self.initUI()

    def initUI(self):
        """Initialize the main user interface settings."""
        self.setWindowTitle('ScanOMetrics')
        self.showMaximized()
        self.setStyleSheet(styles.background_style)

        # Set up the main widget and layout
        self.main_widget = QWidget(self)
        self.setCentralWidget(self.main_widget)
        self.main_layout = QGridLayout(self.main_widget)

        # Add layout separators to the main layout
        self.add_layout_separators()

        # Initialize the main layout areas
        self.title_layout, self.header_layout = QHBoxLayout(), QHBoxLayout()
        self.left_layout, self.right_layout = QVBoxLayout(), QVBoxLayout()

        # Add the layout areas to the main layout grid
        for layout, position in zip([self.title_layout, self.header_layout, self.left_layout, self.right_layout],
                                    [(0, 0), (0, 2), (2, 0), (2, 2)]):
            self.main_layout.addLayout(layout, *position)

        # Set up the title, header, and initial home page
        self.set_title()
        self.set_header()
        self.set_home_page()

    def add_layout_separators(self):
        """Add layout separators (lines) to the main layout for better visual organization."""
        separators_info = [
            (QFrame.VLine, 0, 1, 3, 1),  # Vertical line separator
            (QFrame.HLine, 1, 0, 1, 3)  # Horizontal line separator
        ]

        for frame_shape, row, col, rowspan, colspan in separators_info:
            line = QFrame()
            line.setFrameShape(frame_shape)
            line.setStyleSheet(styles.separator_style)
            self.main_layout.addWidget(line, row, col, rowspan, colspan)

    def set_title(self):
        """Set up the title label at the top of the application."""
        self.title_label = QLabel('Scanometrics')
        self.title_label.setStyleSheet(styles.title_style)
        self.title_label.setTextFormat(Qt.RichText)
        self.title_label.setAlignment(Qt.AlignCenter)
        self.title_layout.addWidget(self.title_label)

    def set_header(self):
        """Set up the header with buttons and other interactive elements."""
        widgets_info = [
            ('info_button', PushButton('Information'), styles.button_style, 'info_page', 'header')
        ]
        self.fill_layout(widgets_info)

    def set_info_page(self):
        """Configure the layout and widgets for the information page."""
        widgets_info = [
            ('left_label', QLabel('Documentation'), "color: #f5f5f5;", None, 'left'),
            ('back_button', PushButton('Back'), styles.button_style, self.current_page, 'left'),
            ('right_label', QLabel('Documentation'), "color: #f5f5f5;", None, 'right')
        ]
        self.fill_layout(widgets_info)

    def set_home_page(self):
        """Configure the layout and widgets for the home page."""
        self.current_page = 'home_page'
        widgets_info = [
            ('left_label', QLabel('Documentation'), "color: #f5f5f5;", None, 'left'),
            ('path_selector', PathSelector(self.sample_folder_path), styles.line_edit_style, 'sample_folder_path',
             'right'),
            ('select_pkl', SelectPKLFile(self.normative_model_string, self.normative_models.keys()), styles.combobox_style, 'normative_model_string', 'right'),
            ('select_sub', SelectSub(self.sub_name), styles.combobox_style, 'sub_name', 'right'),
            ('select_ses', SelectSes(self.ses_nb), styles.combobox_style, 'ses_nb', 'right'),
            ('select_acq_lab', SelectAcqLab(self.acq_lab), styles.combobox_style, 'acq_lab', 'right'),
            ('continue_button', PushButton('right'), styles.button_style, 'continue', 'right')
        ]
        self.fill_layout(widgets_info)

    def set_eval_page(self):
        """Configure the layout and widgets for the evaluation page."""
        self.current_page = 'ev_page'
        widgets_info = [
            ('metric_boxes', MetricCheckBox(), '', 'check_value', 'right'),
            ('select_metrics', SelectROI(self.metric), styles.combobox_style, 'metric', 'left'),
            ('run_proc_button', PushButton('Run Processing Pipeline'), styles.button_style, 'run_proc_pipe', 'left'),
            ('generate_report_button', PushButton('Generate Report'), styles.button_style, 'generate_report', 'left'),
            ('reset_button', PushButton('left'), styles.button_style, 'reset', 'left'),
            ('brain3d', Brain3D(os.path.join(self.sample_folder_path, 'derivatives', self.som_project.metric_proc_pipeline.subjects_dir, "_".join([self.sub_name, self.ses_nb, self.acq_lab]))), '', 'brain_region', 'right')
        ]
        self.fill_layout(widgets_info)
        if self.evaluation_results:
            self.brain3d.output = self.evaluation_results['output']

    def fill_layout(self, widgets_info):
        """Helper method to fill the specified layout with widgets.

        Args:
            widgets_info (list): List of tuples containing widget information:
                                 (widget_name, widget, widget_style, signal_label, layout_position)
        """
        # Mapping layout positions to actual layout objects
        pos_lay = {
            'title': self.title_layout,
            'header': self.header_layout,
            'left': self.left_layout,
            'right': self.right_layout
        }

        # Iterate through widget information and add widgets to the corresponding layout
        for widget_name, widget, widget_style, sgn_label, pos in widgets_info:
            widget.setStyleSheet(widget_style)
            if sgn_label:
                # Connect the widget's signal to the handle_sgn slot, passing the signal label
                widget.sgn.connect(lambda sgn, lbl=sgn_label: self.handle_sgn(sgn, lbl))
            pos_lay[pos].addWidget(widget)
            setattr(self, widget_name, widget)  # Store widget as an attribute of MainWindow

    @pyqtSlot(str)
    def handle_sgn(self, value, signal_type):
        """Handle signals emitted by the widgets and perform corresponding actions."""
        # print(f"Signal received - type: {signal_type}, value: {value}")
        try:
            if signal_type in {'sample_folder_path', 'normative_model_string', 'sub_name', 'ses_nb',
                               'acq_lab', 'check_value', 'metric'}:
                self.get_value(signal_type, value)  # Update the corresponding attribute based on the signal type
            elif signal_type in {'info_page', self.current_page}:
                self.switch_page(signal_type)  # Switch to the requested page
            elif signal_type == 'continue':
                if not self.show_missing_values():  # Check for missing values before proceeding
                    # Set list of subjects to subject of interest, to avoid loading subjects without processed metric values
                    self.som_project.load_subjects(sub_ses_acq_include=[self.som_project.get_subjSesAcq_id(self.sub_name, self.ses_nb, self.acq_lab)])
                    # If data for the selected subject exists, load it and run the evaluation against the normative dataset
                    find_stat_files = os.path.join(self.som_project.metric_proc_pipeline.subjects_dir, "_".join([self.sub_name, self.ses_nb, self.acq_lab]), 'stats2table', '*stats*.txt')
                    found_stat_files = len(glob(find_stat_files))
                    if found_stat_files > 0:
                        self.som_project.load_proc_metrics(ref_metric_values=self.som_project.normativeModel.measured_metrics['orig'],
                                                           ref_metric_names=self.som_project.normativeModel.metric_names,
                                                           ref_covariate_values=self.som_project.normativeModel.covariate_values,
                                                           ref_covariate_names=self.som_project.normativeModel.covariate_names)
                        self.evaluate_subject()
                    self.switch_page('ev_page')
                    self.select_metrics.get_items(self.som_project.metric_names)
            elif signal_type == 'run_proc_pipe':  # Run the processing pipeline
                print('Run Proc Pipeline button pressed')
                print('Running proc pipeline with 2 cpus')
                self.som_project.run_proc_pipeline(subject_id=self.sub_name, n_threads=2)  # TODO: add a way to set n_threads from the GUI
                self.som_project.proc2table(n_threads=1)
                self.som_project.load_proc_metrics(ref_metric_values=self.som_project.normativeModel.measured_metrics['orig'],
                                                   ref_metric_names=self.som_project.normativeModel.metric_names,
                                                   ref_covariate_values=self.som_project.normativeModel.covariate_values,
                                                   ref_covariate_names=self.som_project.normativeModel.covariate_names)
                self.evaluate_subject()
                print("Run Proc Pipeline completed successfully")
            elif signal_type == 'brain_region':  # Show the SpiderPlot of the region depends on metrics (volume, thickness and thicknessstd)
                if value['type'] == 'screenshots':
                    self.brain_screenshots.append(value)
                if value['type'] == 'region_info':
                    self.brain_region_info = value
                    self.figure_dialog = FigureWindow(self)
                    self.figure_dialog.sgn.connect(lambda sgn, lbl='add_fig': self.handle_sgn(sgn, lbl))
            elif signal_type == 'add_fig':
                for fig_data in value:
                    self.plot_data_list.append(fig_data)
            elif signal_type == 'generate_report':
                self.generate_report()
            elif signal_type == 'reset':
                # Confirm and reset the user input variables
                if QMessageBox.question(self, 'Confirm Reset', 'Do you want to reset your entry?',
                                        QMessageBox.Yes | QMessageBox.No, QMessageBox.No) == QMessageBox.Yes:
                    self.reset_variables()
                self.switch_page('home_page')
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}\n\n{traceback.format_exc()}")

    def get_value(self, signal_type, value):
        setattr(self, signal_type, value)
        if signal_type == 'sample_folder_path':
            if self.normative_model_string != "Select your normative model (.pkl)":
                self.initialize_som_project()
                self.som_project.load_normative_model(self.normative_model_string)
            # Update subjects list when sample folder path is selected
            self.select_sub.get_items(self.sample_folder_path)
        elif signal_type == 'normative_model_string':
            # self.normative_model_string = value already set with setattr(self, signal_type, value) at beginning of if loop
            if self.sample_folder_path != "Select the sample folder path...":
                self.initialize_som_project()
                self.som_project.load_normative_model(self.normative_model_string)
        elif signal_type == 'sub_name':
            # Update sessions list when subject is selected
            self.select_ses.get_items(self.sample_folder_path, self.sub_name)
        elif signal_type == 'ses_nb':
            # Update acquisition labels when session is selected
            self.select_acq_lab.get_items(self.sample_folder_path, self.sub_name, self.ses_nb)
        elif signal_type == 'acq_lab':
            self.som_project.load_subjects(sub_ses_acq_include=[self.som_project.get_subjSesAcq_id(self.sub_name, self.ses_nb, self.acq_lab)], covariates_include=['age'])  # Dummy load of subjects without covariates, as default pkl file might not have the proper cov2float dictionary
        elif signal_type == 'check_value':
            self.metric_boxes.reset_checkboxes()  # Uncheck all other checkboxes excepted the actual checked one
            self.brain3d.update_display(self.check_value)  # Update 3D brain display with color
        elif signal_type == 'metric':
            self.figure_dialog = FigureWindow(self, 'cloud_points')
            self.figure_dialog.sgn.connect(lambda sgn, lbl='add_fig': self.handle_sgn(sgn, lbl))

    def switch_page(self, page):
        """Switch to a different page.

        Args:
            page (str): The page to switch to.
        """
        # Remove all widgets from current layout
        self.del_all_widgets()
        # Set the new page based on the input
        {'home_page': self.set_home_page,
         'ev_page': self.set_eval_page,
         'info_page': self.set_info_page}.get(page, lambda: None)()

    def initialize_som_project(self):
        """Initialize the ScanOMetrics project with the selected sample folder and normative model. Here the pkl file
        just has to exist to extract the proc_pipeline_id. Appropriate versioning is taken care of in the call to
        self.som_project.load_normative_model(self.normative_model_string)
        """
        self.som_project = ScanOMetrics_project(self.sample_folder_path, proc_pipeline=self.normative_models[self.normative_model_string]['proc_pipeline_id'], atlas=self.normative_models[self.normative_model_string]['atlas'])

    def evaluate_subject(self):
        """Evaluate the selected subject using the selected metric."""
        try:
            """
            print(f"Evaluating subject: {self.sub_name} for session: {self.ses_nb}"
                  f" with label: {self.acq_lab}")
            """

            # Show the progress dialog
            self.progress_dialog = QProgressDialog("Evaluation in progress ...\n please wait", None, 0, 100, self)
            self.progress_dialog.setWindowModality(Qt.WindowModal)
            self.progress_dialog.setValue(0)
            self.progress_dialog.show()

            # Create and start the evaluation thread
            self.evaluation_thread = EvaluationThread(self.som_project, self.sub_name, self.ses_nb)
            self.evaluation_thread.progress.connect(self.progress_dialog.setValue)
            self.evaluation_thread.finished.connect(self.on_evaluation_finished)
            self.evaluation_thread.error.connect(self.show_error)
            self.evaluation_thread.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error during evaluation: {e}\n\n{traceback.format_exc()}")
            if self.progress_dialog:
                self.progress_dialog.close()

    def on_evaluation_finished(self, subj_covariate_values, normModel_matching_cols, subj_measured_metrics, output,
                               pid):
        """Handle the event when the evaluation is finished."""
        QMessageBox.information(self, "Evaluation Complete", "The evaluation is complete.")
        self.progress_dialog.close()
        # Add the evaluate subject output to the main window and in the Brain3D class
        self.output = self.brain3d.output = output
        """print("subj_covariate_values:", subj_covariate_values)
        print("normModel_matching_cols: ", normModel_matching_cols)
        print("subj_measured_metrics: ", subj_measured_metrics)
        print("Evaluation PID: ", pid)
        print("Evaluation Output: ", self.output)"""

        # Store evaluation results
        self.evaluation_results = {
            "subj_covariate_values": subj_covariate_values,
            "normModel_matching_cols": normModel_matching_cols,
            "subj_measured_metrics": subj_measured_metrics,
            "output": output,
            "pid": pid
        }

    def generate_report(self):
        """
        Generate a PDF report containing the selected plots and metrics.

        This method generates plots for the specified metrics by default and adds them to the report.
        """
        try:
            # Initialize the ReportGenerator
            self.report_generator = ReportGenerator(
                filename=f'{self.sub_name}_report.pdf',
                sub_name=self.sub_name,
                sub_ses=self.ses_nb,
                acq=self.acq_lab,
                folder_path=self.sample_folder_path
            )

            # Extract the necessary values from the evaluation results
            subj_covariate_values = self.evaluation_results.get('subj_covariate_values', [])
            normModel_matching_cols = self.evaluation_results.get('normModel_matching_cols', [])
            subj_measured_metrics = self.evaluation_results.get('subj_measured_metrics', [])
            output = self.evaluation_results.get('output', {})
            pid = self.evaluation_results.get('pid', 0)

            # Define the metrics to be included in the report
            default_metrics = [
                'asegvolume_lhCortexVol',
                'asegvolume_Left-Cerebral-White-Matter',
                'asegvolume_Left-Ventricle-all',
                'asegvolume_Left-Cerebellum',
                'asegvolume_Brain-Stem',
                '%s_lh_MeanThickness_thickness' % aparc_code[self.normative_models[self.normative_model_string]['atlas']]
            ]

            # Generate and collect plots for default metrics
            default_plot_data = []
            for metric in default_metrics:
                fig = self.som_project.plot_single_metric(
                    subject_id=self.sub_name,
                    selected_metric=metric,
                    selected_session=self.ses_nb,
                    selected_acquisition=self.acq_lab,
                    output=output,
                    # subj_covariate_values=subj_covariate_values,
                    # normModel_matching_cols=normModel_matching_cols,
                    # subj_measured_metrics=subj_measured_metrics,
                    matching_covariates=["sex", "sequence"],
                    pid=pid,
                    alpha_uncorr=0.01
                )

                # Save the plot to a BytesIO object
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)

                # Collect the plot data
                default_plot_data.append({'image': buf, 'caption': metric})

            # Add default metrics to the report
            self.report_generator.add_default_metrics(default_plot_data)

            # Add all stored plots to the report generator
            for plot_data in self.plot_data_list:
                self.report_generator.add_plot(plot_data)

            # Add brain screenshots
            for brain_screenshot in self.brain_screenshots:
                self.report_generator.add_brain_snapshot(brain_screenshot)

            # Create the PDF report
            self.report_generator.create_pdf_report()

            # Show a message box indicating the report has been generated
            QMessageBox.information(self, 'Report Generated',
                                    f'Report saved as {self.report_generator.filename}')

            self.is_report_created = True  # Set the flag to True after creating the report

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"An unexpected error occurred while generating the report: {e}\n\n{traceback.format_exc()}")

    def show_error(self, error_msg, traceback_msg):
        """Show error message in case of an exception during the evaluation."""
        QMessageBox.critical(self, "Error", f"Error during evaluation: {error_msg}\n\n{traceback_msg}")
        if self.progress_dialog:
            self.progress_dialog.close()

    def closeEvent(self, event):
        """
        Intercepts the window close event to prompt the user if they attempt to close the window
        without generating a report.

        Parameters:
        event (QCloseEvent): The close event triggered by the user.

        This method will check if a report has been created (is_report_created flag).
        If not, it will display a message box asking for confirmation to close the window.
        If the user chooses 'No', the close event will be ignored.
        If the user chooses 'Yes', the window will be closed.
        """
        if not self.is_report_created:
            # Prompt the user with a message box if the report has not been created
            reply = QMessageBox.question(self, 'Message',
                                         "No report generated, are you sure to exit?",
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.No:
                # Ignore the close event if the user chooses 'No'
                event.ignore()
            else:
                # Accept the close event if the user chooses 'Yes'
                event.accept()
        else:
            # Accept the close event if the report has been created
            event.accept()

    def show_missing_values(self):
        """Display a message box showing which values are missing and return the missing values.

        Returns:
            list: A list of missing values.
        """
        # Check which input values are missing
        missing_values = [
            "Sample folder" if self.sample_folder_path == "Select the sample folder path..." else None,
            "PKL file" if self.normative_model_string == "Select your normative model (.pkl)" else None,
            "Subject" if self.sub_name == "Select the Subject" else None,
            "Session" if self.ses_nb == "Select the Session" else None,
            "Acquisition label" if self.acq_lab == "Select acquisition label" else None
        ]
        missing_values = [value for value in missing_values if value]
        if missing_values:
            # Show warning message if any value is missing
            missing_message = "The following values are missing:\n" + "\n".join(missing_values)
            QMessageBox.warning(self, "Missing Values", missing_message)
        return missing_values

    def del_all_widgets(self):
        """Delete all widgets from both left and right layouts except brain3d."""
        for layout in [self.left_layout, self.right_layout]:
            for i in reversed(range(layout.count())):
                layout.itemAt(i).widget().setParent(None)

    def reset_variables(self):
        """Reset the variables used to get the user entry to their default values."""
        self.sample_folder_path, self.normative_model_string, self.sub_name, self.ses_nb, self.acq_lab, self.metric = \
            ("Select the sample folder path...", "Select your normative model (.pkl)",
             "Select the Subject", "Select the Session", "Select acquisition label", "Select the metrics")
        self.plot_data_list = []  # Reset the plot data list
        self.is_report_created = False


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
