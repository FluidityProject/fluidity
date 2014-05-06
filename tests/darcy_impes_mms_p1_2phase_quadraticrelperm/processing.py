## PARAMETERS 

case_name = "darcy_impes_mms_p1_2phase_quadraticrelperm"
mesh_type_list = ["reg"]                         # [see note 1]
mesh_suffix_list_per_dimension = [
    ["B", "C", "D", "E"],
    ["B", "C", "D"],
    ["B", "C"]]
field_name_list = ["Phase2::Saturation"]         # [2]
norm_list = [2]
mesh_A_number_of_timesteps = 10                  # [2]


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
    case_name, mesh_type_list, mesh_suffix_list_per_dimension, 
    field_name_list, norm_list, mesh_A_number_of_timesteps)

# pass client commands
for arg in argv[1:]:
    test_helper.do(arg)

    
# Notes:
# 
# [1] irreg meshes also possible, but their convergence rates are unruly
# 
# [2] the starting number of timesteps has been increased to 10 because
#     with too few, the rates can look OK even for bad convergence
#     (TODO: check impact on total test time).  However, it has also
#     been discovered that, at certain levels of timestep and mesh
#     refinement, the pressure fields may reconfigure themselves to be
#     completely different to the forced solution.  This destroys the
#     pressure field convergence rates.  Saturation fields seem fine, so
#     only their convergence will be measured.  (TODO: investigate)
