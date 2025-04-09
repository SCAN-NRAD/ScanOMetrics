from pathlib import Path
from reportlab.lib import colors
from reportlab.graphics.shapes import Drawing, Line
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, SimpleDocTemplate, Image, PageBreak, Table, TableStyle
from PIL import Image as pil_image

class ReportGenerator:
    """
    A class to generate a PDF report for ScanOMetrics.

    Attributes:
        filename (str): The name of the PDF file to be generated.
        sub_name (str): The subject name.
        sub_ses (str): The session number.
        acq (str): The acquisition label.
        folder_path (str): The directory to save the report.
        plots (list): A list to store plot data.
        metrics (dict): A dictionary to store metrics data.
        story (list): A list to hold the flowable elements for the PDF.
        styles (StyleSheet1): The default stylesheet for the PDF.
    """

    def __init__(self, filename, sub_name, sub_ses, acq, folder_path):
        """
        Initialize the ReportGenerator with the given parameters.

        Args:
            filename (str): The filename for the generated PDF report.
            sub_name (str): The subject name.
            sub_ses (str): The session number.
            acq (str): The acquisition label.
            folder_path (str): The directory to save the report.
        """
        self.filename = filename
        self.sub_name = sub_name
        self.sub_ses = sub_ses
        self.acq = acq
        self.folder_path = folder_path
        self.plots = []
        self.brain_screenshots = []
        self.metrics = {}
        self.story = []  # List to hold the flowable elements for the PDF
        self.styles = getSampleStyleSheet()  # Get the default style sheet

        # Add a custom caption style
        caption_style = ParagraphStyle(name='Caption', fontSize=10, leading=12)
        self.styles.add(caption_style)

        # Set the filename to the user-specified path
        self.filename = self.get_documents_path()

    def get_documents_path(self):
        """
        Get the path to the user-specified directory and create necessary directories.

        Returns:
            Path: The full path to the PDF file to be generated.
        """
        documents_folder = Path(self.folder_path) / 'ScanOMetrics_report' / self.sub_name / self.sub_ses / self.acq
        # Create the directories if they do not exist
        documents_folder.mkdir(parents=True, exist_ok=True)
        # Return the full path to the PDF file
        return documents_folder / f'{self.sub_name}_{self.sub_ses}_{self.acq}_report.pdf'

    def add_plot(self, plot_data):
        """
        Add a plot to the report.

        Args:
            plot_data (dict): Dictionary containing the plot image and caption.
        """
        self.plots.append({'image': plot_data['image'], 'caption': plot_data.get('caption', '')})

    def add_brain_snapshot(self, brain_screenshot):
        self.brain_screenshots.append(brain_screenshot)

    def add_line(self):
        """
        Add a horizontal line to the report for visual separation.
        """
        self.story.append(Spacer(1, 12))  # Add space after the line
        line_width = 450  # Width of the line (in points)
        line_height = 1  # Height of the line (in points)
        line = Drawing(line_width, line_height)
        line.add(Line(0, 0, line_width, 0, strokeColor=colors.black, strokeWidth=1))
        self.story.append(line)
        self.story.append(Spacer(1, 12))  # Add space after the line

    def add_default_metrics(self, metrics):
        """
        Add the default metrics to the report.

        Args:
            metrics (list): List of default metrics to be added.
        """
        # Add the report title
        self.story.append(Paragraph(f"ScanOMetrics Report for {self.sub_name}", self.styles['Title']))
        self.story.append(Paragraph(f"Session: {self.sub_ses}", self.styles['Heading2']))
        self.story.append(Paragraph(f"Acquisition: {self.acq}", self.styles['Heading2']))
        self.add_line()

        self.story.append(Paragraph("Whole Brain Information", self.styles['Heading2']))
        print("Adding default metric images")
        for metric in metrics:
            if pil_image.open(metric['image']).size[0] == 3000:
                img = Image(metric['image'], width=8 * inch, height=2 * inch)  # Set the image dimensions
            else:
                img = Image(metric['image'], width=8/3 * inch, height=2 * inch)  # Set the image dimensions
            print("%s: %dx%d (originial size was %dx%d)" % (metric['caption'], img._width, img._height, pil_image.open(metric['image']).size[0], pil_image.open(metric['image']).size[1]))
            self.story.append(img)
            self.story.append(Paragraph(metric['caption'], self.styles['Caption']))

    def create_pdf_report(self):
        """
        Create a PDF report with the added plots and metrics.
        """
        # Create a SimpleDocTemplate object with the specified filename and page size
        doc = SimpleDocTemplate(str(self.filename), pagesize=letter)
        styles = self.styles

        self.story.append(PageBreak())
        self.add_line()
        self.story.append(Paragraph("Data Selection", self.styles['Heading2']))

        # Add each plot image and its caption to the PDF
        print("Adding optional plots")
        for plot in self.plots:
            orig_width, orig_height = pil_image.open(plot['image']).size
            img = Image(plot['image'], width=orig_width/orig_height*2 * inch, height=2 * inch)  # Set the image dimensions
            print("%s: %dx%d (originial size was %dx%d)" % (plot['caption'], img._width, img._height, orig_width, orig_height))
            self.story.append(img)  # Add the image to the story
            self.story.append(Paragraph(plot['caption'], styles['Caption']))  # Add the caption to the story

        for brain_screenshot in self.brain_screenshots:
            orig_width, orig_height = pil_image.open(brain_screenshot['left_hemisphere']['image']).size
            colorbar_width, colorbar_height = pil_image.open(brain_screenshot['color_bar']['image']).size
            brain_table = Table([[Image(brain_screenshot['left_hemisphere']['image'], width=orig_width/orig_height*3.5 * inch, height=3.5 * inch),
                                 Image(brain_screenshot['color_bar']['image'], width=colorbar_width/colorbar_height*2*inch, height=2*inch),
                                 Image(brain_screenshot['right_hemisphere']['image'], width=orig_width/orig_height*3.5 * inch, height=3.5 * inch)
                                 ]])
            brain_table.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'MIDDLE')]))
            self.story.append(brain_table)

        # Build the PDF document with the collected elements
        doc.build(self.story)
