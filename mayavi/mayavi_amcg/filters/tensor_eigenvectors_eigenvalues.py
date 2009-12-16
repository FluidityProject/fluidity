## To Do:
## More descriptive description and use to correct terminology

# Author: Daryl Harrison

# Enthought library imports.
from enthought.traits.api import Instance, List, String, Int, Bool
from enthought.traits.ui.api import View, Group, Item

from enthought.tvtk.api import tvtk
from math import *
from numpy import *

# Local imports.
from enthought.mayavi.core.filter import Filter
from enthought.mayavi.core.traits import DEnum

################################################################################
# `TensorEigenvectorsEigenvalues` class.
################################################################################
class TensorEigenvectorsEigenvalues(Filter):
    """
    Computes the eigenvectors and eigenvalues of tensor fields.
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    grid = Instance(tvtk.UnstructuredGrid, allow_none=False)

    processed_tensors_for_eigenvectors_eigenvalues = List()
    processed_tensors_for_edgelengths = List()
    processed_tensors_for_proportionaleigenvectors = List()

    active_tensor = String
    dimensions = Int

    eigenvector_eigenvalue_choice = DEnum(values_name='_dimension_list')
    _dimension_list = List(String)

    edgelength_on = Bool
    proprotionaleigenvectors_on = Bool

    ######################################################################
    # The view.
    ######################################################################
    traits_view = \
        View(
            Group(
                Item(name='eigenvector_eigenvalue_choice', label='Eigenvectors/eigenvalues'),
                Item(name='edgelength_on', label='Edge lengths'),
                Item(name='proprotionaleigenvectors_on', label='Proportional eigenvectors'),
            )
        )

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def update_pipeline(self):
        if len(self.inputs) == 0 or len(self.inputs[0].outputs) == 0:
            return

        self.grid = tvtk.UnstructuredGrid()
        ## This doesn't work - there must be more to copy.
        #self.grid.points = tvtk.Points()
        #self.grid.points.deep_copy(self.inputs[0].outputs[0].points)
        self.grid.deep_copy(self.inputs[0].outputs[0])

        self.dimensions = size(self.grid.points[0])

        for i in range(self.dimensions):
            self._dimension_list.append(i)

        self.active_tensor = self.inputs[0].outputs[0].point_data.tensors.name
        self._set_outputs([self.grid])

    def update_data(self):
        self.active_tensor = self.inputs[0].outputs[0].point_data.tensors.name
        self.data_changed = True

    ######################################################################
    # Non-public interface.
    ######################################################################
    def _eigenvector_eigenvalue_choice_changed(self):
        self._update_output()

    def _edgelength_on_changed(self):
        self._update_output()

    def _proprotionaleigenvectors_on_changed(self):
        self._update_output()

    def _active_tensor_changed(self):
        self.calculate_eigenvectors_eigenvalues()
        self._update_output()

    def _update_output(self):
        if (self.edgelength_on):
            self.calculate_edgelengths()
            self.grid.point_data.set_active_scalars(self.active_tensor+'_edgelengths_'+self.eigenvector_eigenvalue_choice)
        else:
            self.grid.point_data.set_active_scalars(self.active_tensor+'_eigenvalues_'+self.eigenvector_eigenvalue_choice)

        if (self.proprotionaleigenvectors_on):
            self.calculate_propotionaleigenvectors()
            self.grid.point_data.set_active_vectors(self.active_tensor+'_proportionaleigenvectors_'+self.eigenvector_eigenvalue_choice)
        else:
            self.grid.point_data.set_active_vectors(self.active_tensor+'_eigenvectors_'+self.eigenvector_eigenvalue_choice)

        self.pipeline_changed = True

    def calculate_eigenvectors_eigenvalues(self):
        if not (self.active_tensor in self.processed_tensors_for_eigenvectors_eigenvalues):
            self.processed_tensors_for_eigenvectors_eigenvalues.append(self.active_tensor)

            input_grid = self.inputs[0].outputs[0]
            tensor_field = array(input_grid.point_data.get_array(self.active_tensor))

            result = map(lambda x: linalg.eig(reshape(x,(self.dimensions,self.dimensions))), tensor_field)
            sorted_indices = map(lambda x: argsort(x[0]), result)

            for i in range(self.dimensions):
                eigenvalues  = map(lambda x,y: x[0][y[i]], result, sorted_indices)
                eigenvectors = map(lambda x,y: x[1][y[i]], result, sorted_indices)

                eigenvalues_field = tvtk.FloatArray(name=self.active_tensor+'_eigenvalues_'+`i`)
                eigenvalues_field.from_array(eigenvalues)
                self.grid.point_data.add_array(eigenvalues_field)

                eigenvectors_field = tvtk.FloatArray(name=self.active_tensor+'_eigenvectors_'+`i`)
                eigenvectors_field.from_array(eigenvectors)
                self.grid.point_data.add_array(eigenvectors_field)

    def calculate_edgelengths(self):
        if not (self.active_tensor in self.processed_tensors_for_edgelengths):
            self.processed_tensors_for_edgelengths.append(self.active_tensor)

            input_grid = self.inputs[0].outputs[0]

            for i in range(self.dimensions):
                eigenvalues = array(self.grid.point_data.get_array(self.active_tensor+'_eigenvalues_'+`i`))
                edgelengths = map(lambda x: 1/sqrt(x), eigenvalues)

                edgelengths_field = tvtk.FloatArray(name=self.active_tensor+'_edgelengths_'+`i`)
                edgelengths_field.from_array(edgelengths)
                self.grid.point_data.add_array(edgelengths_field)

    def calculate_propotionaleigenvectors(self):
        if not (self.active_tensor in self.processed_tensors_for_proportionaleigenvectors):
            self.processed_tensors_for_proportionaleigenvectors.append(self.active_tensor)

            self.calculate_edgelengths()
            input_grid = self.inputs[0].outputs[0]

            for i in range(self.dimensions):
                eigenvectors = array(self.grid.point_data.get_array(self.active_tensor+'_eigenvectors_'+`i`))
                edgelengths = array(self.grid.point_data.get_array(self.active_tensor+'_edgelengths_'+`i`))
                proportionaleigenvectors = map(multiply, eigenvectors, edgelengths)

                proportionaleigenvectors_field = tvtk.FloatArray(name=self.active_tensor+'_proportionaleigenvectors_'+`i`)
                proportionaleigenvectors_field.from_array(proportionaleigenvectors)
                self.grid.point_data.add_array(proportionaleigenvectors_field)
