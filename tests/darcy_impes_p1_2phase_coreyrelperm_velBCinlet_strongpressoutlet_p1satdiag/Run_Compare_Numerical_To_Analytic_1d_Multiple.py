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

file_name = "darcy_impes_p1_2phase_coreyrelperm_velBCinlet_strongpressoutlet_p1satdiag_relpermupwind_1d"

sol_file_name = [file_name+"_A_17.vtu", 
                 file_name+"_B_33.vtu", 
                 file_name+"_C_65.vtu", 
                 file_name+"_D_129.vtu"]

ana_file_name = "reference_solution/analytic_BL_CoreyPerm_saturation.txt"

out_file_name = [file_name+"_A_17_output", 
                 file_name+"_B_33_output", 
                 file_name+"_C_65_output", 
                 file_name+"_D_129_output"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i])


file_name = "darcy_impes_p1_2phase_coreyrelperm_velBCinlet_strongpressoutlet_p1satdiag_modrelpermupwind_satfesweby_1d"

sol_file_name = [file_name+"_A_17.vtu", 
                 file_name+"_B_33.vtu", 
                 file_name+"_C_65.vtu", 
                 file_name+"_D_129.vtu"]

ana_file_name = "reference_solution/analytic_BL_CoreyPerm_saturation.txt"

out_file_name = [file_name+"_A_17_output", 
                 file_name+"_B_33_output", 
                 file_name+"_C_65_output", 
                 file_name+"_D_129_output"]

for i in range(len(sol_file_name)):
   
   os.system("./Compare_Numerical_To_Analytic_1d.py "+sol_file_name[i]+" "+ana_file_name+" "+out_file_name[i])
