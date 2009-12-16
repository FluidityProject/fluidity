## To Do:
## Test, test, test (in particular, check div - can't tell if output is correct)
## More descriptive description and use correct terminology
## Disable UI when calculating
## Assuming UnstructuredGrid (input and saving) - make this general
## Only works with PointData at present
## Nicer to check ALL preconditions before performing operation, rather than step-by-step
## 'Are you sure?' dialog boxes for clearing form, overwriting field, removing field
## Use of Enum, DEnum for selecting input fields, rather than text input
## Use Collection to add SetActiveAttribute group
##
## Points to Note:
## Grad, Div, Curl, d/dx, d/dy, d/dx all use tvtk.CellDerivatives and tvtk.CellDataToPointData,
## so numerically may not be accurate - also involves a lot of superfluous processing

# Author: Daryl Harrison

# Enthought library imports
from enthought.traits.api import Instance, Enum, Button, List, Tuple, Float, String, File, Code
from enthought.traits.ui.api import View, Group, Item, ListEditor, TupleEditor
from enthought.traits.ui.menu import OKButton

from enthought.tvtk.api import tvtk
from math import *
from numpy import *

# Local imports
from enthought.mayavi.core.filter import Filter
#from enthought.mayavi.core.traits import DEnum

################################################################################
# `FieldOperations` class.
################################################################################
class FieldOperations(Filter):
    """
    Performs standard field operations. Use this filter in conjunction with SetActiveAttribute filter
    to view new/modified fields.
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    grid = Instance(tvtk.UnstructuredGrid, allow_none=False)

    operation       = Enum('Sum','Multiply','Grad','Div','Curl','d/dx','d/dy','d/dz','Normalize','2-norm','Reciprocal','Dot product','Cross product','Step','Natural log','Exponential','Scale/Offset','x-component','y-component','z-component','Python function')
    custom_function = Code('\n# Define a single-line operation in terms of\n# the variable \'inputs\' which is an array of\n# input arrays. math and numpy are imported.\n# The operation must return a single array.\n')

    # Input and output fields
    scale_factor = Float(1.0)
    offset       = Float(0.0)
    field_name   = String

    input_field  = Tuple(field_name, scale_factor, offset)
    input_fields = List(input_field)

    output_field = Tuple(field_name, scale_factor, offset)

    # Buttons
    apply = Button
    clear = Button

    # Saving
    output_file = File
    save        = Button

    # Renaming
    rename_field_old = String
    rename_field_new = String
    rename           = Button

    # Deleting
    remove_field = String
    remove       = Button

    ######################################################################
    # The view.
    ######################################################################

    field_editor      = TupleEditor(labels=['Name','Scale factor','Offset'])
    field_list_editor = ListEditor(style='custom', rows=3, editor=field_editor)

    traits_view = \
        View(
            Group(
                Item(name='operation'),
                Group(
                    Item(name='custom_function'),
                    visible_when='operation==\'Python function\'',
                    show_labels=False
                ),
                Group(
                    Item(name='input_fields', style='custom', editor=field_list_editor),
                    show_labels=False,
                    show_border=True,
                    label='Input fields'
                ),
                Group(
                    Item(name='output_field', style='custom', editor=field_editor),
                    show_labels=False,
                    show_border=True,
                    label='Output field'
                ),
                Group(
                    Item(name='apply', label='Apply operation'),
                    Item(name='clear', label='Clear'),
                    show_labels=False
                ),
                label='Operation'
            ),
            Group(
                Group(
                    Item(name='output_file'),
                    Item(name='save', label='Save'),
                    show_labels=False,
                    show_border=True,
                    label='Save changes to file'
                ),
                Group(
                    Group(
                        Item(name='rename_field_old', label='Field to rename'),
                        Item(name='rename_field_new', label='New name')
                    ),
                    Item(name='rename', label='Rename'),
                    show_labels=False,
                    show_border=True,
                    label='Rename field'
                ),
                Group(
                    Group(
                        Item(name='remove_field', label='Field to remove')
                    ),
                    Item(name='remove', label='Remove'),
                    show_labels=False,
                    show_border=True,
                    label='Remove field'
                ),
                label='Misc'
            ),
            height=600,
            width=550
        )

    dialog_msg = String

    ######################################################################
    # `Filter` interface.
    ######################################################################
    def setup_pipeline(self):
        print 'setup_pipeline'

    def update_pipeline(self):
        print 'update_pipeline'

        if len(self.inputs) == 0 or len(self.inputs[0].outputs) == 0:
            return

        self.grid = tvtk.UnstructuredGrid()
        ## Way of doing this without a deep_copy (idea: copy input to output, with additional arrays tagged on)?
        self.grid.deep_copy(self.inputs[0].outputs[0])

        # Add a blank input field by default
        self.input_fields.append(self.input_field)

        self._set_outputs([self.grid])

    def update_data(self):
        print 'update_data'
        self.data_changed = True

    ######################################################################
    # Non-public interface.
    ######################################################################
    def _apply_fired(self):
        try:
            self.field_operation(self.operation, self.input_fields, self.output_field)
            self._dialog_box('Success', 'Operation complete. The result has been stored in \''+self.output_field[0]+'\'.')
        except Exception, inst:
            self._dialog_box('Error', inst.__str__())

    def _dialog_box(self, title, msg):
        self.dialog_msg = msg
        view = \
            View(
                Group(
                    Item(name='dialog_msg', style='readonly'),
                    show_labels=False
                ),
                title=title,
                buttons=[OKButton],
            )
        self.edit_traits(view,kind='livemodal')

    def _clear_fired(self):
        self.reset_traits(['input_field', 'input_fields', 'output_field'])
        self.input_fields.append(self.input_field)

    def _save_fired(self):
        gridwriter = tvtk.XMLUnstructuredGridWriter()
        gridwriter.file_name = self.output_file
        gridwriter.input = self.grid
        gridwriter.update()

    def _rename_fired(self):
        array = self.grid.point_data.get_array(self.rename_field_old)
        if (array is not None):
            array.name = self.rename_field_new

            self._dialog_box('Success', self.rename_field_old+' has been renamed to '+self.rename_field_new+'.')
            self.reset_traits(['rename_field_old', 'rename_field_new'])
            self.pipeline_changed = True
        else:
            self._dialog_box('Error', self.rename_field_old+' does not exist.')

    def _remove_fired(self):
        if (self.grid.point_data.get_array(self.remove_field) is not None):
            self.grid.point_data.remove_array(self.remove_field)

            self._dialog_box('Success', self.remove_field+' has been removed.')
            self.reset_traits(['remove_field'])
            self.pipeline_changed = True
        else:
            self._dialog_box('Error', self.remove_field+' does not exist.')

    ######################################################################
    # Operations.
    ######################################################################
    def field_operation(self, operation, input_arrays, output_array):
        scalar_op = False       # Operation on scalars only?
        vector_op = False       # Operation on vectors only?
        custom_fn = False       # Custom Python function defined?
        multiple_inputs = True  # Operation with multiple inputs?

        if (operation == 'Sum'):
            self.check_min_input_size(input_arrays, 2)
            function = lambda x: add.reduce(x)
        elif (operation == 'Multiply'):
            self.check_min_input_size(input_arrays, 2)
            function = lambda x: multiply.reduce(x)
        elif (operation == 'Dot product'):
            vector_op = True
            self.check_exact_input_size(input_arrays, 2)
            function = lambda x: map(dot, x[0], x[1])
        elif (operation == 'Cross product'):
            vector_op = True
            self.check_exact_input_size(input_arrays, 2)
            function = lambda x: cross(x[0], x[1])
        elif (operation == 'Python function'):
            custom_fn = True
            try:
                function = lambda inputs: eval(self.custom_function)
            except: 
                raise Exception, 'There is an error in your custom Python function.'

        else:
            multiple_inputs = False
            if (operation == 'Reciprocal'):
                function = lambda x: 1/x
            elif (operation == 'Exponential'):
                function = lambda x: exp(x)
            elif (operation == 'Natural log'):
                function = lambda x: log(x)
            elif (operation == 'Step'):
                function = lambda x: map(lambda x: (x>0)+0, x)
            elif (operation == 'Scale/Offset'):
                function = lambda x: x
            elif (operation == 'd/dx'):
                scalar_op = True
                function = lambda x: self.derivative(x, 0)
            elif (operation == 'd/dy'):
                scalar_op = True
                function = lambda x: self.derivative(x, 1)
            elif (operation == 'd/dz'):
                scalar_op = True
                function = lambda x: self.derivative(x, 2)
            elif (operation == 'Curl'):
                vector_op = True
                function = lambda x: self.grad_curl(x, True)
            elif (operation == 'Grad'):
                scalar_op = True
                function = lambda x: self.grad_curl(x, False)
            elif (operation == 'Div'):
                vector_op = True
                function = self.divergence
            elif (operation == 'Normalize'):
                vector_op = True
                function = lambda x: map(multiply, x, 1/array(map(linalg.norm, x)))
            elif (operation == '2-norm'):
                vector_op = True
                function = lambda x: map(linalg.norm, x)
            elif (operation == 'x-component'):
                vector_op = True
                function = lambda x: self.get_axis(x, 0)
            elif (operation == 'y-component'):
                vector_op = True
                function = lambda x: self.get_axis(x, 1)
            elif (operation == 'z-component'):
                vector_op = True
                function = lambda x: self.get_axis(x, 2)
            else:
                raise Exception, 'Undefined operation.'
            self.check_exact_input_size(input_arrays, 1)

        # Retrieve array data
        input_data = self.retrieve_arrays(input_arrays, scalar_op, vector_op)

        if (multiple_inputs):
            if (custom_fn):
                try:
                    result = function(input_data)
                except:
                    raise Exception, 'There is an error in your custom Python function.'
            else:
                result = function(input_data)
        else:
            result = function(input_data[0])

        output_array_name = output_array[0]
        scale_factor      = output_array[1]
        offset            = output_array[2]

        # Scale and offset result array if necessary
        if (scale_factor != 1.0 or offset != 0.0):
            result = array(result)*scale_factor + offset

        # Store result (replaces array with name output_array_name if it exists)
        output_array = tvtk.FloatArray(name=output_array_name)
        output_array.from_array(result)
        self.grid.point_data.add_array(output_array)
        self.pipeline_changed = True

        ## print array(result)

    def retrieve_arrays(self, input_arrays, scalar_op, vector_op):
        input_data = []
        num_components = None
        for entry in input_arrays:
            array_name   = entry[0]
            scale_factor = entry[1]
            offset       = entry[2]
            arr = self.grid.point_data.get_array(array_name)

            if (arr is None):
                raise Exception, 'Input field \''+array_name+'\' does not exist.'

            if (num_components is None):  # first array we're processing
                num_components = size(arr[0])
                if (vector_op and num_components < 2):
                    raise Exception, 'Input field \''+array_name+'\' is not a vector.'
                if (scalar_op and num_components != 1):
                    raise Exception, 'Input field \''+array_name+'\' is not a scalar.'
            else:
                # check that all input arrays have the same number of components
                if (num_components != size(arr[0])):
                    raise Exception, 'Incompatible input arrays (e.g. scalar and vector).'

            # Scale and offset input array if necessary
            if (scale_factor != 1.0 or offset != 0.0):
                arr = array(arr)*scale_factor + offset

            input_data.append(array(arr))
        return input_data

    def check_min_input_size(self, input_array_names, n):
        if (size(input_array_names,0) < n):
            multiple = ''
            if (n>1):
                multiple = 's'
            raise Exception, 'This operation requires at least '+`n`+' input field'+multiple+'.'

    def check_exact_input_size(self, input_array_names, n):
        if (size(input_array_names,0) != n):
            multiple = ''
            if (n>1):
                multiple = 's'
            raise Exception, 'This operation requires exactly '+`n`+' input field'+multiple+'.'

    def grad_curl(self, array, curl):
        # Can't remove arrays by index - that's why I've given it a name
        temp_name  = 'a_temp_name_hopefully_nobody_will_ever_use'
        temp_array = tvtk.FloatArray(name=temp_name)
        temp_array.from_array(array)
        self.grid.point_data.add_array(temp_array)

        if(curl):
            self.grid.point_data.set_active_vectors(temp_name)
        else:
            self.grid.point_data.set_active_scalars(temp_name)

        cd = tvtk.CellDerivatives()
        cd.input = self.grid
        if (curl):
            cd.vector_mode = 'compute_vorticity'
        cd.update()

        cp = tvtk.CellDataToPointData()
        cp.input = cd.output
        cp.update()

        self.grid.point_data.remove_array('a_temp_name_hopefully_nobody_will_ever_use')

        return cp.output.point_data.get_array('Vorticity')

    def divergence(self, array):
        x_component = self.get_axis(array, 0)
        y_component = self.get_axis(array, 1)
        z_component = self.get_axis(array, 2)

        derivative_x = self.derivative(x_component, 0)
        derivative_y = self.derivative(y_component, 1)
        derivative_z = self.derivative(z_component, 2)

        return add.reduce([derivative_x, derivative_y, derivative_z])

    def derivative(self, array, axis):
        array = self.grad_curl(array, False)
        return self.get_axis(array, axis)

    def get_axis(self, array, axis):
        return array[:,axis]