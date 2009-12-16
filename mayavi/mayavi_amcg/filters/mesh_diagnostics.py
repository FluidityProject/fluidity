## To Do:
## Make stats fields read-only
## More descriptive description

# Author: Daryl Harrison

# Enthought library imports
from enthought.traits.api import Instance, Enum, Bool, Tuple, Float
from enthought.traits.ui.api import View, Item, Group, TupleEditor
from enthought.tvtk.api import tvtk

# Local imports
from enthought.mayavi.core.filter import Filter
from enthought.mayavi.core.traits import DEnum
from enthought.mayavi.core.dataset_manager import DatasetManager
from numpy import array
from sets import Set

################################################################################
# `MeshDiagnostics` class.
################################################################################
class MeshDiagnostics(Filter):
    """
    Identifies stiff elements of a solid mesh.
    Also includes diagnostics from vtkMeshQuality:
    <http://www.vtk.org/doc/release/5.0/html/a01739.html>
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    _dataset_manager = Instance(DatasetManager, allow_none=False)
    _surface_grid = Instance(tvtk.PointSet, allow_none=False)

    stiff_elements = Bool
    _first_time_stiff_elements = Bool(True)

    _mesh_quality_filter = Instance(tvtk.MeshQuality, args=(), allow_none=False)

    tet_quality_measure = Enum('Edge ratio', 'Aspect ratio', 'Radius ratio', 'Min angle', 'Frobenius norm', 'Volume')
    triangle_quality_measure = Enum('Edge ratio', 'Aspect ratio', 'Radius ratio', 'Min angle', 'Frobenius norm')
    quad_quality_measure = Enum('Edge ratio', 'Aspect ratio', 'Radius ratio', 'Min angle', 'Med Frobenius norm', 'Max Frobenius norm')
    hex_quality_measure = Enum('Edge ratio')

    tet_stats = triangle_stats = quad_stats = hex_stats = Tuple(Float,Float,Float,Float,Float)

    ######################################################################
    # The view.
    ######################################################################
    stats_editor = TupleEditor(labels=['Minimum','Average','Maximum','Variance','Number of cells'])
    traits_view = \
        View(
            Item(name='stiff_elements'),
            Group(
                Group(
                    Item(name='tet_quality_measure'),
                    Item(name='tet_stats', style='custom', editor=stats_editor),
                    label='Tetrahedron quality measure',
                    show_labels=False,
                    show_border=True,
                ),
                Group(
                    Item(name='triangle_quality_measure'),
                    Item(name='triangle_stats', style='custom', editor=stats_editor),
                    label='Triangle quality measure',
                    show_labels=False,
                    show_border=True,
                ),
                Group(
                    Item(name='quad_quality_measure'),
                    Item(name='quad_stats', style='custom', editor=stats_editor),
                    label='Quadrilateral quality measure',
                    show_labels=False,
                    show_border=True,
                ),
                Group(
                    Item(name='hex_quality_measure'),
                    Item(name='hex_stats', style='custom', editor=stats_editor),
                    label='Hexahedron quality measure',
                    show_labels=False,
                    show_border=True,
                ),
                enabled_when='stiff_elements==False',
                show_border=False
            ),
            width=300
        )

    """
    Hexahedron quality measure    - EdgeRatio

    Quadrilateral quality measure - EdgeRatio
                                  - AspectRatio
                                  - RadiusRatio
                                  - MinAngle
                                  - MedFrobeniusNorm
                                  - MaxFrobeniusNorm

    Tetrahedron quality measure   - EdgeRatio
                                  - AspectRatio
                                  - RadiusRatio
                                  - MinAngle
                                  - FrobeniusNorm
                                  - Volume (compatibility mode on/off)

    Triangle quality measure      - EdgeRatio
                                  - AspectRatio
                                  - RadiusRatio
                                  - MinAngle
                                  - FrobeniusNorm
    """

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def update_pipeline(self):
        if len(self.inputs) == 0 or len(self.inputs[0].outputs) == 0:
            return

        self._mesh_quality_filter.set_input(self.inputs[0].outputs[0])
        for type in ['tet','triangle','quad','hex']:
            self._change_quality_measure(type, 'Edge ratio')
        self._mesh_quality_filter.update()
        self._set_outputs([self._mesh_quality_filter.output])

    def update_data(self):
        self.data_changed = True

    ######################################################################
    # Non-public interface.
    ######################################################################

    def _stiff_elements_changed(self, value):
        if (value):
            self._find_stiff_elements()
            self._set_outputs([self._surface_grid])
        else:
            self._set_outputs([self._mesh_quality_filter.output])

    def _tet_quality_measure_changed(self, value):
        self._change_quality_measure('tet', value)

    def _change_quality_measure(self, cell_type, quality_measure):
        if (quality_measure == 'Volume'):
            self._mesh_quality_filter.compatibility_mode = True
            self._mesh_quality_filter.volume = True
        else:
            self._mesh_quality_filter.compatibility_mode = False
            self._mesh_quality_filter.volume = False
            quality_measure_names = {'Edge ratio':'edge_ratio', 'Aspect ratio':'aspect_ratio', 'Radius ratio':'radius_ratio',
                   'Min angle':'min_angle', 'Frobenius norm':'frobenius_norm', 'Max Frobenius norm':'max_frobenius_norm',
                   'Med Frobenius norm':'med_frobenius_norm'}
            setattr(self._mesh_quality_filter, '%s_quality_measure' %cell_type, quality_measure_names.get(quality_measure))

        self._mesh_quality_filter.update()
        self._set_outputs([self._mesh_quality_filter.output])

        cell_type_fullnames = {'tet':'Tetrahedron', 'triangle':'Triangle', 'quad':'Quadrilateral', 'hex':'Hexahedron'}
        stats = self._mesh_quality_filter.output.field_data.get_array('Mesh %s Quality' %cell_type_fullnames.get(cell_type))[0]
        setattr(self, '%s_stats' %cell_type, stats)

    def _find_stiff_elements(self):
        if (self._first_time_stiff_elements):
            self._first_time_stiff_elements = False
            self._surface_grid = type(self.inputs[0].outputs[0])()
            self._surface_grid.copy_structure(self.inputs[0].outputs[0])

            if (self._surface_grid.points.data_type == 'double'):
                # tvtk.DataSetSurfaceFilter produces surface points as float
                # so we need to convert input points from doubles to floats
                # if necessary
                float_points = array(self._surface_grid.points, 'float32')
                self._surface_grid.points = float_points

            self._dataset_manager = DatasetManager(dataset=self._surface_grid)
            surface_points = self._get_surface_points()
            self._add_stiff_elements(surface_points)

    def _get_surface_points(self):
        surface_filter = tvtk.DataSetSurfaceFilter()
        surface_filter.set_input(self._surface_grid)
        surface_filter.update()
        return surface_filter.output.points

    def _add_stiff_elements(self, surface_points):
        stiff_elements = []
        surface_points = Set(surface_points)

        for i in range(self._surface_grid.number_of_cells):
            cell = self._surface_grid.get_cell(i)
            points = Set(cell.points)
            val = int(points.issubset(surface_points))
            stiff_elements.append(val)

        array_name = 'Stiff Elements'
        self._dataset_manager.add_array(array(stiff_elements), array_name, 'cell')
        self._dataset_manager.activate(array_name, 'cell')