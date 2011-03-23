SHELL=sh

input: clean
	
clean:
	rm -f gmon.out  input.dat  mirror_int_data.dat  time_count
	$(MAKE) -C BL_test_10elements_1d clean
	$(MAKE) -C BL_test_10elements_1d_DG clean
	$(MAKE) -C Sodds_shock_tube_1phase_1d clean
	$(MAKE) -C BL_test_gravity clean
	$(MAKE) -C BL_test_10elements_1d_compressible clean
