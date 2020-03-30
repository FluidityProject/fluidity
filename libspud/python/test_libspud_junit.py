import os
import libspud

from junit_xml import TestSuite, TestCase

suite = TestSuite('libspud for python',[])

dirpath = os.path.dirname(os.path.abspath(__file__))

def test(test_str, exception=None):
    """ Pass/fail test."""
    suite.test_cases.append(TestCase('libspud for python.%s'%test_str))    
    try:
        result = eval(test_str)
        if not result:
             suite.test_cases[-1].add_failure_info('Failure')
    except Exception as e:
        suite.test_cases[-1].add_failure_info('Exception', str(e))

def exception_test(test_str, exception):
    """Test should throw exception."""
    try:
        suite.test_cases.append(TestCase("%s fails"%test_str))
        eval(test_str)
        suite.test_cases[-1].add_failure_info('No exception')
    except exception as e:
        return
    # reach here on test failure
    suite.test_cases[-1].add_failure_info('Exception', str(e))

try:
    suite.test_cases.append(TestCase('libspud for python.load_options'))
    libspud.load_options(dirpath+'/test.flml')
except:
    suite.test_cases[-1].add_failure_info('Exception')
    with open('test_results.xml', 'w') as handle:
        suite.to_file(handle, [suite])
    raise ValueError

test("libspud.get_option('/timestepping/timestep') == 0.025")
test("libspud.get_number_of_children('/geometry') == 5")
test("libspud.get_child_name('geometry', 0) == 'dimension'")

test("libspud.option_count('/problem_type') == 1")
test("libspud.have_option('/problem_type')")

test("libspud.get_option_type('/geometry/dimension') is int")
test("libspud.get_option_type('/problem_type') is str")

test("libspud.get_option_rank('/geometry/dimension') == 0")
test("libspud.get_option_rank('/physical_parameters/gravity/vector_field::GravityDirection/prescribed/value/constant') == 1")

test("libspud.get_option_shape('/geometry/dimension') == (-1, -1)")
test("libspud.get_option_shape('/problem_type')[0] > 1")
test("libspud.get_option_shape('/problem_type')[1] == -1")

test("libspud.get_option('/problem_type') == 'multimaterial'")
test("libspud.get_option('/geometry/dimension') == 2")
libspud.set_option('/geometry/dimension', 3)


test("libspud.get_option('/geometry/dimension') == 3")

list_path = '/material_phase::Material1/scalar_field::MaterialVolumeFraction/prognostic/boundary_conditions::LetNoOneLeave/surface_ids'
test("libspud.get_option_shape(list_path) == (4, -1)")
test("libspud.get_option_rank(list_path) == 1")
test("libspud.get_option(list_path) == [7, 8, 9, 10]")

libspud.set_option(list_path, [11, 12, 13, 14, 15])
test("libspud.get_option_shape(list_path) == (5, -1)")
test("libspud.get_option_rank(list_path) == 1")
test("libspud.get_option(list_path)==[11, 12, 13, 14, 15]")

tensor_path = '/material_phase::Material1/tensor_field::DummyTensor/prescribed/value::WholeMesh/anisotropic_asymmetric/constant'
test("libspud.get_option_shape(tensor_path) == (2, 2)")
test("libspud.get_option_rank(tensor_path) == 2")

test("libspud.get_option(tensor_path)==[[1.0,2.0],[3.0,4.0]]")

libspud.set_option(tensor_path, [[5.0,6.0,2.0],[7.0,8.0,1.0]])
test("libspud.get_option_shape(tensor_path) == (2,3)")
test("libspud.get_option_rank(tensor_path) == 2")

test("libspud.get_option(tensor_path)==[[5.0, 6.0, 2.0],[7.0, 8.0, 1.0]]")

exception_test("libspud.add_option('/foo')", libspud.SpudNewKeyWarning)

test("libspud.option_count('/foo') == 1")

libspud.set_option('/problem_type', 'helloworld')
test("libspud.get_option('/problem_type') == 'helloworld'")

exception_test("libspud.set_option_attribute('/foo/bar', 'foobar')", libspud.SpudNewKeyWarning)

test("libspud.get_option('/foo/bar') == 'foobar'")
  
libspud.delete_option('/foo')
test("libspud.option_count('/foo') == 0")

exception_test("libspud.get_option('/foo')", libspud.SpudKeyError)
exception_test("libspud.get_option('/geometry')", libspud.SpudTypeError)

libspud.write_options(dirpath+'test_out.flml')

exception_test("libspud.set_option('/test', 4.3)", libspud.SpudNewKeyWarning)
test("libspud.get_option('/test') == 4.3")
test("libspud.set_option('/test',4) or True")

test("libspud.get_option('/test') == 4")

libspud.set_option('/test',[[4.0,2.0,3.0],[2.0,5.0,6.6]]) 

test("libspud.get_option('/test') == [[4.0,2.0,3.0],[2.0,5.0,6.6]]")

libspud.set_option('/test',"Hallo")

test("libspud.get_option('/test') == 'Hallo'")

libspud.set_option('/test',[1,2,3])

test("libspud.get_option('/test') == [1,2,3]")

libspud.set_option('/test',[2.3,3.3])

test("libspud.get_option('/test') == [2.3,3.3]")

exception_test("libspud.set_option('/test')", libspud.SpudError)

  
with open(os.path.dirname(os.path.abspath(__file__))+'/test_results.xml', 'w') as handle:
    suite.to_file(handle, [suite])
