"""This filter enables one to select a portion of an input dataset
using a plane and clip it so only one side remains.

Many thanks to Prabhu for ScalarCutPlane and TransformData.
Wouldn't have been able to code this otherwise.

"""
# Author: Samir Talwar <samir.talwar06@imperial.ac.uk>
# License: BSD Style.

# Enthought library imports.
from enthought.traits.api import Instance, Int, Trait, TraitMap, Button
from enthought.traits.ui.api import View, Group, Item

from enthought.tvtk.api import tvtk

# Local imports
from enthought.mayavi.core.common import error
from enthought.mayavi.filters.filter_base import FilterBase
from enthought.mayavi.core.pipeline_info import PipelineInfo
from enthought.mayavi.components.implicit_plane import ImplicitPlane


######################################################################
# `DataSetClipper` class.
######################################################################
class DataSetClipper(FilterBase):
    # The version of this class.  Used for persistence.
    __version__ = 0
    
    # The implicit plane widget used to easily place the implicit function.
    implicit_plane = Instance(ImplicitPlane, allow_none=False)
    
    # The actual filter.
    filter = Instance(tvtk.ClipDataSet, allow_none=False)
    
    # I'm not sure what this'll work with. vtkUnstructuredGrid is confirmed.
    # Everything else is somewhat possible.
    input_info = PipelineInfo(datasets=['any'],
                              attribute_types=['any'],
                              attributes=['any'])
    
    output_info = PipelineInfo(datasets=['any'],
                               attribute_types=['any'],
                               attributes=['any'])
    
    ########################################
    # View related traits.
    
    # Button to reset the boundaries of the plane.
    # This should really be done automatically.
    reset_button = Button('Reset Boundaries')
    
    # The View for this object.
    view = View(Group(Item(name='reset_button'),
                      Item(name='implicit_plane', style='custom'),
                      show_labels=False,
                      ),
                )
    
    ######################################################################
    # `Filter` interface
    ######################################################################
    def setup_pipeline(self):
        self.implicit_plane = ImplicitPlane()
        self.filter = tvtk.ClipDataSet()
    
    def update_pipeline(self):
        inputs = self.inputs
        if len(inputs) == 0:
            return
        
        implicit_plane = self.implicit_plane
        implicit_plane.inputs = inputs
        implicit_plane.update_pipeline()
        
        widget = self.implicit_plane.widget
        widget.outline_translation = 0
        self.widgets = [widget]
        
        filter = self.filter
        filter.input = inputs[0].outputs[0]
        filter.clip_function = implicit_plane.plane
        filter.update()
        
        self._set_outputs([filter.output])
        
        self.pipeline_changed = True
    
    def update_data(self):
        # Do nothing if there is no input.
        if len(self.inputs) == 0:
            return
        
        # Propagate the data_changed event.
        self.data_changed = True
    
    ######################################################################
    # Non-public methods.
    ######################################################################
    def _on_implicit_plane_changed(self):
        self.filter.clip_function = self.implicit_plane.plane
        self.filter.update()
        self.render()
    
    def _reset_button_fired(self):
        if len(self.widgets) == 0:
            return
        
        self.widgets[0].place_widget()
        self.render()

