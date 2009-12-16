## To Do:
## More descriptive description and use to correct terminology
##
## Points to Note:
## Filter only works if it lies directly below the source in the pipeline. It only operates on point data.

# Original code by Tim Bond <http://amcg.ese.ic.ac.uk/index.php?title=Local:Using_Mayavi2>
# Author: Daryl Harrison

# Enthought library imports.
from enthought.traits.api import Instance, Float, Array, List, String
from enthought.traits.ui.api import View, Group, Item
from enthought.tvtk.api import tvtk

# Local imports
from enthought.mayavi.core.common import debug
from enthought.mayavi.core.filter import Filter
from numpy import linalg, array, sum, sqrt, column_stack, arctan2, arccos, zeros, cross, vdot


######################################################################
# `ProjectionAndDepthStretch` class.
######################################################################
class ProjectionAndDepthStretch(Filter):
    """
    Depth stretching and simple projection of earth's surface onto a plane.
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    scale_factor = Float(0.0001)
    grid = Instance(tvtk.UnstructuredGrid, allow_none=False)

    longitude = Array
    colatitude = Array
    initial_depth = Array
    unit_r  = Array
    unit_lon = Array
    unit_lat = Array
    processed_vectors = List

    active_scalar = String
    active_vector = String

    # The view
    view = \
        View(
            Group(
                Item(name='scale_factor', style='simple'),
            ),
            scrollable=True,
            resizable=True,
        )

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def setup_pipeline(self):
        debug('setup_pipeline')
        """Override this method so that it *creates* its tvtk
        pipeline.

        This method is invoked when the object is initialized via
        `__init__`.  Note that at the time this method is called, the
        tvtk data pipeline will *not* yet be setup.  So upstream data
        will not be available.  The idea is that you simply create the
        basic objects and setup those parts of the pipeline not
        dependent on upstream sources and filters.
        """
        # Scene is updated when scale_factor trait changed.
        #self.scale_factor.on_trait_change(self.render)

    def update_pipeline(self):         # Called when we drag filter under another VTK file
        debug('update_pipeline')
        """Override this method so that it *updates* the tvtk pipeline
        when data upstream is known to have changed.

        This method is invoked (automatically) when the input fires a
        `pipeline_changed` event.
        """
        # Do nothing if there is no input.
        if len(self.inputs) == 0:
            return

        magn = linalg.norm
        earth_radius = 6378000.0

        # By default we set the input to the first output of the first input
        self.grid = tvtk.UnstructuredGrid()
        self.grid.deep_copy(self.inputs[0].reader.output)  ## WAY OF DOING THIS WITHOUT A DEEP COPY?
        #self.inputs[0].outputs[0] ## DOESN'T WORK WITH update_data()

        # Split array by column into the cartesian coordinates
        xyz = array(self.grid.points)
        x = xyz[:,0]
        y = xyz[:,1]
        z = xyz[:,2]

        xyz_squared = xyz**2
        r_squared = xyz_squared.sum(axis=1)
        r = sqrt(r_squared)

        self.longitude = arctan2(y, x)
        self.colatitude = 1.3 * arccos(z/r)
        self.initial_depth = earth_radius - r

        # Now we do vector correction
        self.unit_r = column_stack((x/r, y/r, z/r))
        len_lon = sqrt(x**2 + y**2)
        self.unit_lon = column_stack((y/len_lon, -1*x/len_lon, zeros(len(x))))
        self.unit_lat = cross(self.unit_r, self.unit_lon)

        # Correct current vector array
        current_vector_array = self.inputs[0].outputs[0].point_data.vectors.name
        self._correct_vector_array(current_vector_array)

        self._apply_stretch()

        # Propagate the data_changed event - let modules that depend on this filter know that pipeline has changed
        self.pipeline_changed = True

    def update_data(self):             # Called when we change an option under VTK file, e.g. one of the vectors
        debug('update_data')
        """Override this method to do what is necessary when upstream
        data changes.

        This method is invoked (automatically) when any of the inputs
        sends a `data_changed` event.
        """
        # Do nothing if there is no input.
        if len(self.inputs) == 0:
            return

        if (self.inputs[0].reader.output.point_data.scalars):
            self.active_scalar = self.inputs[0].reader.output.point_data.scalars.name
        if (self.inputs[0].reader.output.point_data.vectors):
            self.active_vector = self.inputs[0].reader.output.point_data.vectors.name

        # Propagate the data_changed event - let modules that depend on this filter know that data has changed
        self.data_changed = True

    def _active_scalar_changed(self, value):
        self.grid.point_data.set_active_scalars(value)
        self._update_display()

    def _active_vector_changed(self, value):
        self.grid.point_data.set_active_vectors(value)
        self._update_display()

    def _scale_factor_changed(self, old, new):
        debug('_scale_factor_changed')
        self._apply_stretch()

    def _apply_stretch(self):
        depth = self.initial_depth * self.scale_factor
        self.grid.points = column_stack((self.longitude, self.colatitude, depth))
        self._update_display()

    def _update_display(self):
        self._set_outputs([self.grid])

    def _correct_vector_array(self, vector_array):
        if not (vector_array in self.processed_vectors):
            self.processed_vectors.append(vector_array)

            vector = self.grid.point_data.get_array(vector_array) # Reference to vector array
            xyz_vector = array(vector)

            mag_lon = array(map(vdot, xyz_vector, self.unit_lon))
            mag_lat = map(vdot, xyz_vector, self.unit_lat)
            mag_r = map(vdot, xyz_vector, self.unit_r)

            vector.from_array(column_stack((-1*mag_lon, mag_lat, mag_r)))