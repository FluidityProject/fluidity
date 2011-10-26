SHELL=sh

input: clean
	
clean:
	rm -f gmon.out  input.dat  mirror_int_data.dat  time_count
	$(MAKE) -C Advection_DomeHat_1D clean
	$(MAKE) -C Advection_PointedHat_1D clean
	$(MAKE) -C Advection_TopHat_1D clean
	$(MAKE) -C BL_1d clean
	$(MAKE) -C BL_1d_Compressible clean
	$(MAKE) -C BL_1d_DGPress clean
	$(MAKE) -C BL_1d_Gravity clean
	$(MAKE) -C BL_1d_GravityUnstable clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_LVI_NonEquilibrium_DG clean
	$(MAKE) -C CompModel_1d_2Phases_2Components_NonEquilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_2Components clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_Equilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_NonEquilibrium clean
	$(MAKE) -C CompModel_1d_2Phases_3Components_LVI_NonEquilibrium clean
	$(MAKE) -C shocktube_1d clean

