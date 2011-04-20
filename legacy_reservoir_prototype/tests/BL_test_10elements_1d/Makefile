SHELL = sh

SIM = BL_test_1d

input: clean

	../../../scripts/interval --dx=0.2 0 4 1d_mesh

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
	gfortran -o input_satura_BL1dg_func  input_satura_BL1dg_func.f90
	gfortran -o input_perm_BL1_func      input_perm_BL1_func.f90
	gfortran -o input_eos_coef_BL1_func  input_eos_coef_BL1_func.f90
	gfortran -o input_usource_BL1_func   input_usource_BL1_func.f90

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
	rm -f input_satura_BL1dg_func        
	rm -f input_perm_BL1_func        
	rm -f input_eos_coef_BL1_func       
	rm -f input_usource_BL1_func
	rm -f gmon.out
	rm -f *.d.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
