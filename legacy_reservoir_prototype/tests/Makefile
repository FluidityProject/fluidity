SHELL=sh

input: clean
	
clean:
	rm -f gmon.out  input.dat  mirror_int_data.dat  time_count
	$(MAKE) -C Sodds_shock_tube_1phase_1d clean
	$(MAKE) -C Advection_1D clean
	$(MAKE) -C BL_test_10elements_1d_compressible clean
	$(MAKE) -C BL_1d clean
	$(MAKE) -C BL_1d_DGPress clean
	$(MAKE) -C BL_1d_Gravity clean
	$(MAKE) -C BL_1d_GravityUnstable clean
	$(MAKE) -C BL_test_10elements_2d clean 
	$(MAKE) -C CompModel_1d_2Phases_3Components_LVI_NonEquilibrium_DG clean
	$(MAKE) -C CompModel_1d_2Phases_2Components_NonEquilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_2Components clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_Equilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_NonEquilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_LVI_NonEquilibrium clean

