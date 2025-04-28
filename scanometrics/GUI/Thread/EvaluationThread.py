import traceback
from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np


class EvaluationThread(QThread):
    """
    A QThread subclass for evaluating a subject.

    This thread handles the process of evaluating a subject's session using the provided
    metric. It emits signals to update the progress, signal completion, and report errors.

    Attributes:
        progress (pyqtSignal): Signal to emit progress updates.
        finished (pyqtSignal): Signal to emit when the evaluation is finished.
        error (pyqtSignal): Signal to emit when an error occurs.
    """
    # Signals to communicate with the main thread
    progress = pyqtSignal(int)
    finished = pyqtSignal(np.ndarray, np.ndarray, dict, dict, float)
    error = pyqtSignal(str, str)

    def __init__(self, som_project, sub_name, ses_nb):
        """
        Initializes the EvaluationThread.

        Args:
            som_project (ScanOMetrics_project): The ScanOMetrics project instance.
            sub_name (str): The name of the subject to evaluate.
            ses_nb (str): The session number of the subject.
        """
        super().__init__()
        self.som_project = som_project
        self.sub_name = sub_name
        self.ses_nb = ses_nb

    def run(self):
        """
        Runs the evaluation process in a separate thread.

        This method simulates the evaluation process by emitting progress signals at
        different stages and finally emits the finished signal with the evaluation results.
        """
        try:
            # Simulate evaluation process with delays and progress updates
            self.sleep(1)  # Simulate some processing time
            self.progress.emit(25)  # Emit progress update (25%)

            # Simulate evaluating the subject
            subj_covariate_values, normModel_matching_cols, subj_measured_metrics, output, pid = self.som_project.evaluate_singleSubject_allSes(
                self.sub_name,
                matching_covariates=['sex', 'sequence'],  # TODO: add option to set covariates in GUI
            )

            self.progress.emit(50)  # Emit progress update (50%)

            self.sleep(1)  # Simulate some more processing time
            self.progress.emit(75)  # Emit progress update (75%)

            self.sleep(1)  # Simulate final processing time
            self.progress.emit(100)  # Emit progress update (100%)

            self.finished.emit(subj_covariate_values, normModel_matching_cols, subj_measured_metrics, output, pid)  # Emit the finished signal with the results

        except Exception as e:
            # Handle exceptions and emit an error signal
            self.error.emit(str(e), traceback.format_exc())  # Emit the error signal with the error message and traceback
