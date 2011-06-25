SHELL = sh

SIM = Sodds_shock_tube_1phase_1d

input: clean

	../../../scripts/interval --dx=0.1 0 1 1d_mesh

	@echo "compute input for $(SIM)"
	gfortran -o input_wic_vol_bc_BL1_fcn input_wic_vol_bc_BL1_fcn.f90
	gfortran -o input_wic_d_bc_BL1_fcn   input_wic_d_bc_BL1_fcn.f90
	gfortran -o input_wic_u_bc_BL1_fcn   input_wic_u_bc_BL1_fcn.f90
	gfortran -o input_wic_p_bc_BL1_fcn   input_wic_p_bc_BL1_fcn.f90

	gfortran -o input_suf_vol_bc_BL1_fcn input_suf_vol_bc_BL1_fcn.f90
	gfortran -o input_suf_d_bc_BL1_fcn   input_suf_d_bc_BL1_fcn.f90
	gfortran -o input_suf_t_bc_BL1_fcn   input_suf_t_bc_BL1_fcn.f90
	gfortran -o input_suf_u_bc_BL1_fcn   input_suf_u_bc_BL1_fcn.f90

	gfortran -o input_satura_BL1_func    input_satura_BL1_func.f90
	gfortran -o input_perm_BL1_func      input_perm_BL1_func.f90
	gfortran -o input_eos_coef_BL1_func  input_eos_coef_BL1_func.f90

	gfortran -o pressure_in              pressure_in.f90

clean:
	rm -f input_wic_vol_bc_BL1_fcn
	rm -f input_wic_d_bc_BL1_fcn       
	rm -f input_wic_u_bc_BL1_fcn        
	rm -f input_wic_p_bc_BL1_fcn        
	rm -f input_suf_vol_bc_BL1_fcn        
	rm -f input_suf_d_bc_BL1_fcn        
	rm -f input_suf_t_bc_BL1_fcn        
	rm -f input_suf_u_bc_BL1_fcn       
	rm -f input_satura_BL1_func        
	rm -f input_perm_BL1_func        
	rm -f input_eos_coef_BL1_func
	rm -f pressure_in 
	rm -f gmon.out
	rm -f *.d.* fluidity.* *.stat 
	rm -f time_count tempfile CheckSparsityMatrix.dat fvalues mirror* fort.*
	rm -f *.ele *.bound *.face *.node
