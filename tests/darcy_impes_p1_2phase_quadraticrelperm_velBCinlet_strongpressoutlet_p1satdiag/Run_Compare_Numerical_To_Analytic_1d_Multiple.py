#!/usr/bin/env python

""" 

 Run_Compare_Numerical_To_Analytic_1d_Multiple.py
 
 This small script will run the script
 Compare_Numerical_To_Analytic_1d.py 
 for a set of hard coded input/ouput 
 filenames associated with this test case.
 
"""

import Compare_Numerical_To_Analytic_1d
import os

# Compare the phase2 saturation for relpermupwind
file_name = "darcy_impes_p1_2phase_quadraticrelperm_velBCinlet_strongpressoutlet_p1satdiag_relpermupwind_1d"

sol_file_name = [file_name+"_A_1.vtu", 
                 file_name+"_B_1.vtu", 
                 file_name+"_C_1.vtu", 
                 file_name+"_D_1.vtu"]

ana_file_name = "reference_solution/analytic_BL_QuadraticPerm_saturation.txt"

out_file_name = [file_name+"_A_1_output_saturation", 
                 file_name+"_B_1_output_saturation", 
                 file_name+"_C_1_output_saturation", 
                 file_name+"_D_1_output_saturation"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i]+" Phase2::Saturation")


# Compare the phase2 saturation for modrelpermupwind_satfesweby
file_name = "darcy_impes_p1_2phase_quadraticrelperm_velBCinlet_strongpressoutlet_p1satdiag_modrelpermupwind_satfesweby_1d"

sol_file_name = [file_name+"_A_1.vtu", 
                 file_name+"_B_1.vtu", 
                 file_name+"_C_1.vtu", 
                 file_name+"_D_1.vtu"]

ana_file_name = "reference_solution/analytic_BL_QuadraticPerm_saturation.txt"

out_file_name = [file_name+"_A_1_output_saturation", 
                 file_name+"_B_1_output_saturation", 
                 file_name+"_C_1_output_saturation", 
                 file_name+"_D_1_output_saturation"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i]+" Phase2::Saturation")

# Compare the phase1 pressure for relpermupwind
file_name = "darcy_impes_p1_2phase_quadraticrelperm_velBCinlet_strongpressoutlet_p1satdiag_relpermupwind_1d"

sol_file_name = [file_name+"_A_1.vtu", 
                 file_name+"_B_1.vtu", 
                 file_name+"_C_1.vtu", 
                 file_name+"_D_1.vtu"]

ana_file_name = "reference_solution/analytic_BL_QuadraticPerm_pressure.txt"

out_file_name = [file_name+"_A_1_output_pressure", 
                 file_name+"_B_1_output_pressure", 
                 file_name+"_C_1_output_pressure", 
                 file_name+"_D_1_output_pressure"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i]+" Phase1::Pressure")


# Compare the phase1 pressure for modrelpermupwind_satfesweby
file_name = "darcy_impes_p1_2phase_quadraticrelperm_velBCinlet_strongpressoutlet_p1satdiag_modrelpermupwind_satfesweby_1d"

sol_file_name = [file_name+"_A_1.vtu", 
                 file_name+"_B_1.vtu", 
                 file_name+"_C_1.vtu", 
                 file_name+"_D_1.vtu"]

ana_file_name = "reference_solution/analytic_BL_QuadraticPerm_pressure.txt"

out_file_name = [file_name+"_A_1_output_pressure", 
                 file_name+"_B_1_output_pressure", 
                 file_name+"_C_1_output_pressure", 
                 file_name+"_D_1_output_pressure"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i]+" Phase1::Pressure")
