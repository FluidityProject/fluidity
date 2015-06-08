from sys import path
path.append('../darcy_impes_common')

from buckley_leverett_test_tools import BuckleyLeverettTestSuite

numerical_filename_stem = "darcy_impes_p1_2phase_quadraticrelperm_velBCinlet_strongpressoutlet_withgravity_updip"
model_name_list = ["relpermupwind", "modrelpermupwind_satfesweby"]
grid_name_list_per_dimension = [["A", "B", "C", "D"], 
                                ["A", "B", "C"],
                                ["A", "B"]]
field_name_list = ["Phase2::Saturation"]
analytic_filename_stem = "reference_solution/analytic_BL_QuadraticPerm_withgravity_updip"

test_suite = BuckleyLeverettTestSuite(
    numerical_filename_stem, model_name_list, grid_name_list_per_dimension,
    analytic_filename_stem, field_name_list)

test_suite.process_folder()
test_suite.generate_reports()
