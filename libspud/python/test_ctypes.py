import libspud

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

try:
  libspud.add_option('/foo')
except libspud.SpudNewKeyWarning:
  pass
print libspud.option_count('/foo')

libspud.set_option('/problem_type', 'helloworld')
print libspud.get_option('/problem_type')

try:
  libspud.set_option_attribute('/foo/bar', 'foobar')
except libspud.SpudNewKeyWarning:
  pass
print libspud.get_option('/foo/bar')
  
libspud.delete_option('/foo')
print libspud.option_count('/foo')

libspud.write_options('test_out.flml')
