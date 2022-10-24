# To Do:
# More descriptive description and use to correct terminology
# Author: Daryl Harrison

import numpy as np
from traits.api import Bool, Instance, Int, List, String
from traitsui.api import Group, Item, View
from tvtk.api import tvtk

from mayavi.core.filter import Filter
from mayavi.core.trait_defs import DEnum

# Local imports.
# Enthought library imports.


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

    eigenvector_eigenvalue_choice = DEnum(values_name="_dimension_list")
    _dimension_list = List(String)

    edgelength_on = Bool
    proprotionaleigenvectors_on = Bool

    ######################################################################
    # The view.
    ######################################################################
    traits_view = View(
        Group(
            Item(
                name="eigenvector_eigenvalue_choice", label="Eigenvectors/eigenvalues"
            ),
            Item(name="edgelength_on", label="Edge lengths"),
            Item(name="proprotionaleigenvectors_on", label="Proportional eigenvectors"),
        )
    )

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def update_pipeline(self):
        if len(self.inputs) == 0 or len(self.inputs[0].outputs) == 0:
            return

        self.grid = tvtk.UnstructuredGrid()
        # This doesn't work - there must be more to copy.
        # self.grid.points = tvtk.Points()
        # self.grid.points.deep_copy(self.inputs[0].outputs[0].points)
        self.grid.deep_copy(self.inputs[0].outputs[0])

        self.dimensions = self.grid.points[0].size

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
        if self.edgelength_on:
            self.calculate_edgelengths()
            self.grid.point_data.set_active_scalars(
                self.active_tensor
                + "_edgelengths_"
                + self.eigenvector_eigenvalue_choice
            )
        else:
            self.grid.point_data.set_active_scalars(
                self.active_tensor
                + "_eigenvalues_"
                + self.eigenvector_eigenvalue_choice
            )

        if self.proprotionaleigenvectors_on:
            self.calculate_propotionaleigenvectors()
            self.grid.point_data.set_active_vectors(
                self.active_tensor
                + "_proportionaleigenvectors_"
                + self.eigenvector_eigenvalue_choice
            )
        else:
            self.grid.point_data.set_active_vectors(
                self.active_tensor
                + "_eigenvectors_"
                + self.eigenvector_eigenvalue_choice
            )

        self.pipeline_changed = True

    def calculate_eigenvectors_eigenvalues(self):
        if not (
            self.active_tensor in self.processed_tensors_for_eigenvectors_eigenvalues
        ):
            self.processed_tensors_for_eigenvectors_eigenvalues.append(
                self.active_tensor
            )

            input_grid = self.inputs[0].outputs[0]
            tensor_field = np.array(input_grid.point_data.get_array(self.active_tensor))

            result = [
                np.linalg.eig(np.reshape(x, (self.dimensions, self.dimensions)))
                for x in tensor_field
            ]
            sorted_indices = [np.argsort(x[0]) for x in result]

            for i in range(self.dimensions):
                eigenvalues = list(map(lambda x, y: x[0][y[i]], result, sorted_indices))
                eigenvectors = list(
                    map(lambda x, y: x[1][y[i]], result, sorted_indices)
                )

                eigenvalues_field = tvtk.FloatArray(
                    name=self.active_tensor + "_eigenvalues_" + repr(i)
                )
                eigenvalues_field.from_array(eigenvalues)
                self.grid.point_data.add_array(eigenvalues_field)

                eigenvectors_field = tvtk.FloatArray(
                    name=self.active_tensor + "_eigenvectors_" + repr(i)
                )
                eigenvectors_field.from_array(eigenvectors)
                self.grid.point_data.add_array(eigenvectors_field)

    def calculate_edgelengths(self):
        if not (self.active_tensor in self.processed_tensors_for_edgelengths):
            self.processed_tensors_for_edgelengths.append(self.active_tensor)

            for i in range(self.dimensions):
                eigenvalues = np.array(
                    self.grid.point_data.get_array(
                        self.active_tensor + "_eigenvalues_" + repr(i)
                    )
                )
                edgelengths = [1 / np.sqrt(x) for x in eigenvalues]

                edgelengths_field = tvtk.FloatArray(
                    name=self.active_tensor + "_edgelengths_" + repr(i)
                )
                edgelengths_field.from_array(edgelengths)
                self.grid.point_data.add_array(edgelengths_field)

    def calculate_propotionaleigenvectors(self):
        if (
            self.active_tensor
            not in self.processed_tensors_for_proportionaleigenvectors
        ):
            self.processed_tensors_for_proportionaleigenvectors.append(
                self.active_tensor
            )

            self.calculate_edgelengths()

            for i in range(self.dimensions):
                eigenvectors = np.array(
                    self.grid.point_data.get_array(
                        self.active_tensor + "_eigenvectors_" + repr(i)
                    )
                )
                edgelengths = np.array(
                    self.grid.point_data.get_array(
                        self.active_tensor + "_edgelengths_" + repr(i)
                    )
                )
                proportionaleigenvectors = list(
                    map(np.multiply, eigenvectors, edgelengths)
                )

                proportionaleigenvectors_field = tvtk.FloatArray(
                    name=self.active_tensor + "_proportionaleigenvectors_" + repr(i)
                )
                proportionaleigenvectors_field.from_array(proportionaleigenvectors)
                self.grid.point_data.add_array(proportionaleigenvectors_field)
