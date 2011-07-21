SHELL = sh

SIM = Advection_1D

input: clean

	../../../scripts/interval --dx=0.05 0 1 1d_mesh

	@echo "compute input for $(SIM)"
	gfortran -o input_suf_t_bc           input_suf_t_bc.f90
	gfortran -o input_suf_one_bc         input_suf_one_bc.f90
	gfortran -o input_wic_t_bc           input_wic_t_bc.f90

clean:
	rm -f input_suf_t_bc  
	rm -f input_wic_t_bc
	rm -f input_suf_one_bc 
	rm -f gmon.out
	rm -f *.d.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat fluidity.*
	rm -f *.ele *.bound *.face *.node
