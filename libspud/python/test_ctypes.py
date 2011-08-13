import libspud
print libspud.__file__

libspud.load_options('test.flml')

libspud.print_options()

print libspud.number_of_children('/geometry')
print libspud.get_child_name('geometry', 0)

print libspud.option_count('/problem_type')
print libspud.have_option('/problem_type')

print libspud.get_option_type('/geometry/dimension')
print libspud.get_option_type('/problem_type')

print libspud.get_option_rank('/geometry/dimension')
print libspud.get_option_rank('/physical_parameters/gravity/vector_field::GravityDirection/prescribed/value/constant')

print libspud.get_option_shape('/geometry/dimension')
print libspud.get_option_shape('/problem_type')

print libspud.get_option('/problem_type')
print libspud.get_option('/geometry/dimension')
libspud.set_option('/geometry/dimension', 3)
print libspud.get_option('/geometry/dimension')

list_path = '/material_phase::Material1/scalar_field::MaterialVolumeFraction/prognostic/boundary_conditions::LetNoOneLeave/surface_ids'
print libspud.get_option_shape(list_path)
print libspud.get_option_rank(list_path)
print libspud.get_option(list_path)
assert(libspud.get_option(list_path)==[7,8,9,10])
libspud.set_option(list_path, [11,12,13,14])
print libspud.get_option_shape(list_path)
print libspud.get_option_rank(list_path)
print libspud.get_option(list_path)
assert(libspud.get_option(list_path)==[11,12,13,14])

tensor_path = '/material_phase::Material1/tensor_field::DummyTensor/prescribed/value::WholeMesh/anisotropic_asymmetric/constant'
print libspud.get_option_shape(tensor_path)
print libspud.get_option_rank(tensor_path)
print libspud.get_option(tensor_path)
assert(libspud.get_option(tensor_path)==[[1.0,2.0],[3.0,4.0]])
libspud.set_option(tensor_path, [[5.0,6.0],[7.0,8.0]])
print libspud.get_option_shape(tensor_path)
print libspud.get_option_rank(tensor_path)
print libspud.get_option(tensor_path)
assert(libspud.get_option(tensor_path)==[[5.0,6.0],[7.0,8.0]])

try:
  libspud.add_option('/foo')
except libspud.SpudNewKeyWarning, e:
  print "caught libspud.SpudNewKeyWarning: "+e.message
print libspud.option_count('/foo')

libspud.set_option('/problem_type', 'helloworld')
print libspud.get_option('/problem_type')

try:
  libspud.set_option_attribute('/foo/bar', 'foobar')
except libspud.SpudNewKeyWarning, e:
  print "caught libspud.SpudNewKeyWarning: "+e.message
print libspud.get_option('/foo/bar')
  
libspud.delete_option('/foo')
print libspud.option_count('/foo')

try:
  libspud.get_option('/foo')
except libspud.SpudKeyError, e:
  print "caught libspud.SpudKeyError: "+e.message

try:
  libspud.get_option('/geometry')
except libspud.SpudTypeError, e:
  print "caught libspud.SpudTypeError: "+e.message

libspud.write_options('test_out.flml')
