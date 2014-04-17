## PARAMETERS 

case_name = "darcy_impes_mms_p1_2phase_quadraticrelperm"
mesh_type_list = ["reg", "irreg"]
mesh_suffix_list_per_dimension = [
    ["A", "B", "C", "D", "E"],
    ["A", "B", "C", "D"],
    ["A", "B", "C"]]
field_name_list = ["Phase1::Pressure", "Phase2::Saturation"]
norm_list = [2]
finish_time = 1.0


## SCRIPT

from sys import argv, exit, path
path.append('../darcy_impes_common')

from manufactured_solution_test_tools import ManufacturedSolutionTestSuite

if len(argv)==1:
    print("""\nRun this script with one or more of the following arguments:
    gen - generate source terms etc. using solution_generator.py
    pre - preprocess (generate meshes and options files)
    proc - process (run simulations)
    post - postprocess (write convergence rates)
    all - all of the above""")
    # TODO put in 'clean' argument
    exit()

# regenerate symbols?
arg1 = str.lower(argv[1][0:3])
if arg1=='gen' or arg1=='all':
    from solution_generator import generate
    generate()

# initialise helper object
from solution_expressions import solution_dict, domain_extents
test_helper = ManufacturedSolutionTestSuite(
    case_name, solution_dict, finish_time, domain_extents,
    mesh_type_list, mesh_suffix_list_per_dimension, field_name_list, norm_list)

# pass client commands
for arg in argv[1:]:
    test_helper.do(arg)
