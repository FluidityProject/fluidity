SHELL=sh

input: clean
	
clean:
	rm -f gmon.out  input.dat  mirror_int_data.dat  time_count
	$(MAKE) -C Advection_SinglePhase clean
	$(MAKE) -C Advection_SinglePhase_Adaptive clean
	$(MAKE) -C Buckley_Leverett_2D clean
	$(MAKE) -C Compositional_Model clean
	$(MAKE) -C Sod_Shock_Tube clean
	$(MAKE) -C Sloshing_Water_Tank clean
	$(MAKE) -C Sod_Shock_Tube_2D clean

