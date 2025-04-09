import io

import nibabel.freesurfer.io as fs
from PyQt5 import QtWidgets, QtCore
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np


class Hemisphere(QtWidgets.QWidget):
    sgn = QtCore.pyqtSignal(dict)
    camera_updated = QtCore.pyqtSignal(vtk.vtkCamera)

    def __init__(self, brain_model_path, annot_file_path, init_camera):
        """Initialize the class with the paths to the brain model and annotation files."""
        super().__init__()
        self.brain_model_path = brain_model_path
        self.annot_file_path = annot_file_path
        self.initial_camera_position = init_camera[0]
        self.initial_camera_focal_point = init_camera[1]
        self.initial_camera_view_up = init_camera[2]

        # Initialize the UI and VTK variables
        self.layout = self.vtkWidget = self.renderer = \
            self.interactor = self.actor = self.region_labels = \
            self.region_colors = self.region_names = self.polydata = None

        # Set up the user interface
        self.setup_ui()
        # Set up the VTK environment
        self.setup_vtk()
        # Load the brain model
        self.load_brain_model()
        # Load the annotation file
        self.load_annot_file()
        # Initialize the interactor
        self.initialize_interactor()

    def setup_ui(self):
        """Set up the user interface layout."""
        self.layout = QtWidgets.QVBoxLayout(self)  # Use a vertical box layout for the widget

        self.vtkWidget = QVTKRenderWindowInteractor(self)  # Create a VTK render window interactor widget
        self.vtkWidget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Expanding)  # Allow the widget to expand

        self.layout.addWidget(self.vtkWidget)  # Add the VTK widget to the layout

    def setup_vtk(self):
        """Set up the VTK renderer and interactor."""
        self.renderer = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtkWidget.GetRenderWindow().GetInteractor()

        interactor_style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(interactor_style)
        self.interactor.AddObserver("EndInteractionEvent", self.on_camera_update)

        # Reset view
        self.reset_camera_position()

    def reset_camera_position(self):
        """Reset the camera to its initial position."""
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(self.initial_camera_position)
        camera.SetFocalPoint(self.initial_camera_focal_point)
        camera.SetViewUp(self.initial_camera_view_up)
        self.renderer.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()

    def on_camera_update(self, caller, event):
        camera = self.renderer.GetActiveCamera()
        self.camera_updated.emit(camera)

    def update_camera(self, camera):
        current_camera = self.renderer.GetActiveCamera()

        # Mirror the camera position
        mirrored_position = [-camera.GetPosition()[0], camera.GetPosition()[1], camera.GetPosition()[2]]
        current_camera.SetPosition(mirrored_position)

        # Mirror the camera focal point
        mirrored_focal_point = [-camera.GetFocalPoint()[0], camera.GetFocalPoint()[1], camera.GetFocalPoint()[2]]
        current_camera.SetFocalPoint(mirrored_focal_point)

        current_camera.SetViewUp(camera.GetViewUp())

        self.renderer.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()

    def load_brain_model(self):
        """Load the brain model from the file and add it to the renderer."""
        vertices, faces = fs.read_geometry(self.brain_model_path)  # Read the geometry from the brain model file

        points = vtk.vtkPoints()  # Create a VTK points object
        for vertex in vertices:  # Add each vertex to the points object
            points.InsertNextPoint(vertex[0], vertex[1], vertex[2])

        polygons = vtk.vtkCellArray()  # Create a VTK cell array for the polygons
        for face in faces:  # Add each face to the cell array
            polygons.InsertNextCell(3)  # Each face is a triangle (3 points)
            polygons.InsertCellPoint(face[0])
            polygons.InsertCellPoint(face[1])
            polygons.InsertCellPoint(face[2])

        polydata = vtk.vtkPolyData()  # Create a VTK polydata object
        polydata.SetPoints(points)  # Set the points of the polydata
        polydata.SetPolys(polygons)  # Set the polygons of the polydata

        self.polydata = polydata  # Store the polydata object

        mapper = vtk.vtkPolyDataMapper()  # Create a VTK polydata mapper
        mapper.SetInputData(self.polydata)  # Set the polydata as input for the mapper

        self.actor = vtk.vtkActor()  # Create a VTK actor
        self.actor.SetMapper(mapper)  # Set the mapper for the actor

        self.actor.GetProperty().SetColor(0.5, 0.5, 0.5)  # RGB values for gray color

        self.renderer.AddActor(self.actor)  # Add the actor to the renderer
        self.renderer.ResetCamera()  # Reset the camera to fit the model in the view

    def load_annot_file(self):
        """Load the annotation file for brain region labeling and coloring."""
        labels, ctab, names = fs.read_annot(self.annot_file_path)  # Read the annotation file
        self.region_labels = labels  # Store the region labels
        self.region_colors = ctab[:, :4]  # Store the region colors
        self.region_names = names  # Store the region names

    def initialize_interactor(self):
        """Initialize the interactor and add event observers."""
        self.interactor.Initialize()  # Initialize the interactor
        self.interactor.Start()  # Start the interactor
        self.interactor.AddObserver("LeftButtonPressEvent", self.click)  # Add an observer for left button clicks

    def click(self, obj, event):
        """Handle mouse click events to select parts of the brain model."""
        click_pos = self.interactor.GetEventPosition()  # Get the position of the click
        picker = vtk.vtkCellPicker()  # Create a cell picker
        picker.Pick(click_pos[0], click_pos[1], 0, self.renderer)  # Perform the pick operation
        clicked_actor = picker.GetActor()  # Get the actor that was clicked

        if clicked_actor == self.actor:  # If the clicked actor is the brain model
            picked_position = picker.GetPickPosition()  # Get the position of the picked point
            point_id = picker.GetPointId()  # Get the ID of the picked point
            if point_id != -1:  # If a valid point was picked
                label = self.region_labels[point_id]  # Get the label of the region
                region_name = self.region_names[label].decode('utf-8')  # Get the name of the region
                region_info = {'type': 'region_info', 'region_name': region_name, 'picked_position': picked_position}
                self.sgn.emit(region_info)  # Emit the signal with region information
            else:
                region_info = f"Selected Area: {picked_position}, Region: Unknown"
                self.sgn.emit(region_info)  # Emit the signal with region information

    def update_polydata_with_annot(self, devtn_logps):
        """Update the polydata with colors from the annotation file."""
        # Create a VTK array to store colors
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(4)  # We need RGBA (4 components)
        colors.SetName("Colors")  # Name the array

        # Initialize all points with a default color (gray)
        num_points = self.polydata.GetNumberOfPoints()
        for _ in range(num_points):
            colors.InsertNextTuple4(80, 80, 80, 255)  # Default gray color

        # Update colors based on region indices and logp values
        for region_idx, region_name in enumerate(self.region_names):
            region_points = np.where(self.region_labels == region_idx)[0]  # Get points belonging to this region
            logpval = devtn_logps[region_idx]
            color = self.logpval_to_color(logpval)
            if color is not None:
                for point in region_points:
                    colors.SetTuple4(point, color[0], color[1], color[2], color[3])

        self.polydata.GetPointData().SetScalars(colors)  # Set the colors for the polydata points
        self.vtkWidget.GetRenderWindow().Render()  # Render the updated polydata

    @staticmethod
    def logpval_to_color(logpval):
        """
        Map log-p-values to colors, creating a gradient similar to the provided image.
        Blue for negative values, gray for values near zero, and red to yellow for positive values.
        """
        # Set to gray if |logpval| < 1.3010299956639813
        plot_thr = 1.3010299956639813
        if -plot_thr < logpval < plot_thr:
            return [80, 80, 80, 255]
        if logpval <= -plot_thr:
            # Blue to Light Blue gradient for negative values less than -1.3
            # Normalize to [0, 1] where -2.0 maps to 1 and -1.3 maps to 0
            normalized_logpval = np.min([-logpval - plot_thr, 2])/2.0
            blue = int(255)
            green = int(normalized_logpval * 255)
            return [0, green, blue, 255]
        if logpval >= plot_thr:
            # Red to Yellow gradient for positive values greater than 1.3
            # Normalize to [0, 1] where 1.3 maps to 0 and 2.0 maps to 1
            normalized_logpval = np.min([logpval - plot_thr, 2])/2.0
            red = int(255)
            green = int(normalized_logpval * 255)
            return [red, green, 0, 255]

    def save_screenshot(self):
        """
        Capture the VTK render window as an image and return it as a dictionary with the image buffer and caption.
        """
        window_to_image_filter = vtk.vtkWindowToImageFilter()
        window_to_image_filter.SetInput(self.vtkWidget.GetRenderWindow())
        window_to_image_filter.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(window_to_image_filter.GetOutputPort())

        # writer.WriteToMemoryOn()
        writer.SetWriteToMemory(True)
        writer.Write()

        data = writer.GetResult()
        buffer = io.BytesIO()
        buffer.write(memoryview(data))
        buffer.seek(0)

        return {'image': buffer}
