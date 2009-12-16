## To Do:
## If a value outside of the original range is applied, the scalar bar range is not updated.
## Refers to TriangleWriter in _save_fired() - change to correct location
## More descriptive description

# Author: Daryl Harrison

# Enthought library imports
from enthought.traits.api import Instance, Range, Bool, Float, List, String, File, Button
from enthought.traits.ui.api import View, Group, Item, ListEditor

from enthought.tvtk.api import tvtk
from enthought.mayavi.core.dataset_manager import DatasetManager
from numpy import *

# Local imports
from enthought.mayavi.core.filter import Filter

################################################################################
# `BoundaryMarkerEditor` class.
################################################################################
class BoundaryMarkerEditor(Filter):
    """
    Edit the boundary marker of a Triangle surface mesh. To use: select the label to
    assign, hover your cursor over the cell you wish to edit, and press 'p'.
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    _current_grid = Instance(tvtk.UnstructuredGrid, allow_none=False)
    _input_grid = Instance(tvtk.UnstructuredGrid, args=(), allow_none=False)
    _extract_cells_filter = Instance(tvtk.ExtractCells, args=(), allow_none=False)
    _dataset_manager = Instance(DatasetManager, allow_none=False)
    _cell_mappings = List

    label_to_apply = Range(0,255)
    select_coplanar_cells = Bool
    epsilon = Range(0.0, 1.0, 0.0001)
    mask_labels = Bool
    labels_to_mask = List(label_to_apply)

    # Saving file
    output_file = File
    save = Button

    ######################################################################
    # The view.
    ######################################################################
    traits_view = \
        View(
            Group(
                Item(name='label_to_apply'),
                Item(name='select_coplanar_cells'),
                Item(name='epsilon', enabled_when='select_coplanar_cells', label='Tolerance'),
                Item(name='mask_labels'),
                Group(
                    Item(name='labels_to_mask', style='custom', editor=ListEditor(rows=3)),
                    show_labels=False,
                    show_border=True,
                    label='Labels to mask',
                    enabled_when='mask_labels==True'
                ),
                Group(
                    Item(name='output_file'),
                    Item(name='save', label='Save'),
                    show_labels=False,
                    show_border=True,
                    label='Save changes to file (give only a basename, without the file extension)'
                )
            ),
            height=500,
            width=600
        )

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def update_pipeline(self):
        if len(self.inputs) == 0 or len(self.inputs[0].outputs) == 0:
            return

        # Call cell_picked() when a cell is clicked.
        self.scene.picker.cellpicker.add_observer("EndPickEvent", self.cell_picked)
        self.scene.picker.pick_type = 'cell_picker'
        self.scene.picker.tolerance = 0.0
        self.scene.picker.show_gui = False

        self._input_grid.deep_copy(self.inputs[0].outputs[0])

        self._current_grid = self._input_grid
        self._dataset_manager = DatasetManager(dataset=self._input_grid)
        self._set_outputs([self._current_grid])

        # Filter for masking.
        self._extract_cells_filter.set_input(self._input_grid)

    def update_data(self):
        self.data_changed = True

    ######################################################################
    # Non-public interface.
    ######################################################################

    def cell_picked(self, object, event):
        cell_id = self.scene.picker.cellpicker.cell_id
        self.modify_cell(cell_id, self.label_to_apply)

        if (self.select_coplanar_cells):
            self.modify_neighbouring_cells(cell_id)

        if (self.mask_labels):
            self.perform_mask()

        self._dataset_manager.activate(self._input_grid.cell_data.scalars.name, 'cell')
        self._dataset_manager.update()
        self.pipeline_changed = True

    def get_all_cell_neigbours(self, cell_id, cell):
        neighbour_cell_ids = array([], dtype=int)

        for i in range(cell.number_of_edges):
            # Get points belonging to ith edge
            edge_point_ids = cell.get_edge(i).point_ids

            # Find neigbours which share the edge
            current_neighbour_cell_ids = tvtk.IdList()
            self._current_grid.get_cell_neighbors(cell_id, edge_point_ids, current_neighbour_cell_ids)
            neighbour_cell_ids = append(neighbour_cell_ids, array(current_neighbour_cell_ids))

        return neighbour_cell_ids.tolist()

    def modify_neighbouring_cells(self, cell_id):
        cell = self._current_grid.get_cell(cell_id)

        cell_normal = [0,0,0]
        cell.compute_normal(cell.points[0], cell.points[1], cell.points[2], cell_normal)

        cells_pending = self.get_all_cell_neigbours(cell_id, cell)
        cells_visited = [cell_id]

        while (len(cells_pending) > 0):
            current_cell_id = cells_pending.pop()

            if (current_cell_id not in cells_visited):
                cells_visited.append(current_cell_id)
                current_cell = self._current_grid.get_cell(current_cell_id)

                current_cell_normal = [0,0,0]
                current_cell.compute_normal(current_cell.points[0], current_cell.points[1], current_cell.points[2], current_cell_normal)

                if (dot(cell_normal, current_cell_normal) > (1-self.epsilon)):
                    self.modify_cell(current_cell_id, self.label_to_apply)
                    cells_pending.extend(self.get_all_cell_neigbours(current_cell_id, current_cell))

    def _mask_labels_changed(self):
        if (self.mask_labels):
            self.perform_mask()
            self._current_grid = self._extract_cells_filter.output
        else:
            self._current_grid = self._input_grid

        self._set_outputs([self._current_grid])
        self.pipeline_changed = True

    def _labels_to_mask_changed(self):
        self.perform_mask()

    def _labels_to_mask_items_changed(self):
        self.perform_mask()

    def perform_mask(self):
        labels_array = self._input_grid.cell_data.get_array(self._input_grid.cell_data.scalars.name)
        in_masked = map(lambda x: x in self.labels_to_mask, labels_array)

        unmasked_cells_list = tvtk.IdList()
        cell_ids = range(self._input_grid.number_of_cells)
        # _cell_mappings is indexed by cell_id of the original input grid, and each value
        # is the new cell_id of the corresponding cell in the masked grid
        self._cell_mappings = map(lambda masked,cell_id: None if masked else unmasked_cells_list.insert_next_id(cell_id), in_masked, cell_ids)

        self._extract_cells_filter.set_cell_list(unmasked_cells_list)
        self._extract_cells_filter.update()
        self.pipeline_changed = True

    def modify_cell(self, cell_id, value):
        if (self.mask_labels):
            cell_id = self._cell_mappings.index(cell_id) # Adjust cell_id if masked
        self._input_grid.cell_data.get_array(self._input_grid.cell_data.scalars.name)[cell_id] = value


    def _save_fired(self):
        from mayavi_amcg.triangle_writer import TriangleWriter
        if (self.output_file):
            writer = TriangleWriter(self._input_grid, self.output_file)
            writer.write()
            print "#### Saved ####"