## PARAMETERS 

case_name = "darcy_impes_mms_p1_2phase_quadraticrelperm"
# irreg meshes also possible, but their convergence rates are unruly
mesh_type_list = ["reg"]
mesh_suffix_list_per_dimension = [
    ["B", "C", "D", "E"],
    ["B", "C", "D"],
    ["B", "C"]]
field_name_list = ["Phase1::Pressure", "Phase2::Saturation"]
norm_list = [2]
finish_time = 1.0


## SCRIPT

from sys import argv, exit, path
path.append('../darcy_impes_common')

if len(argv)==1:
    print("""\nRun this script with one or more of the following arguments:
    xml - make an xml file for the test harness
    gen - generate source terms etc. using solution_generator.py
    pre - preprocess (generate meshes and options files)
    proc - process (run simulations)
    post - postprocess (write convergence rates)
    all - all of the above
    clean - clean everything no longer needed by the test harness
            (*.geo, *.msh, *.diml, *.vtu)""")
    exit()


# initialise helper object
from manufactured_solution_test_tools import ManufacturedSolutionTestSuite
test_helper = ManufacturedSolutionTestSuite(
    case_name, finish_time, mesh_type_list, mesh_suffix_list_per_dimension, 
    field_name_list, norm_list)

# pass client commands
for arg in argv[1:]:
    test_helper.do(arg)
