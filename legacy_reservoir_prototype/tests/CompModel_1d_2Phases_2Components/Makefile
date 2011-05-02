SHELL = sh

SIM = InputComp2P2C_

input: clean

	../../../scripts/interval --dx=0.1 0 1 1d_mesh

	@echo "compute input for $(SIM)"
	gfortran -o $(SIM)KComp   $(SIM)KComp.f90
	gfortran -o $(SIM)permability   $(SIM)permability.f90
	gfortran -o $(SIM)saturation    $(SIM)saturation.f90
	gfortran -o $(SIM)suf_t_bc      $(SIM)suf_t_bc.f90
	gfortran -o $(SIM)suf_vol_bc    $(SIM)suf_vol_bc.f90
	gfortran -o $(SIM)wic_d_bc      $(SIM)wic_d_bc.f90
	gfortran -o $(SIM)wic_p_bc      $(SIM)wic_p_bc.f90
	gfortran -o $(SIM)wic_vol_bc    $(SIM)wic_vol_bc.f90
	gfortran -o $(SIM)EOSCoeff      $(SIM)EOSCoeff.f90
	gfortran -o $(SIM)Comp          $(SIM)Comp.f90
	gfortran -o $(SIM)Comp2		$(SIM)Comp2.f90
	gfortran -o $(SIM)suf_comp_bc   $(SIM)suf_comp_bc.f90

clean:
	rm -f $(SIM)KComp 
	rm -f $(SIM)permability
	rm -f $(SIM)saturation 
	rm -f $(SIM)suf_t_bc
	rm -f $(SIM)suf_vol_bc 
	rm -f $(SIM)wic_d_bc
	rm -f $(SIM)wic_p_bc
	rm -f $(SIM)wic_vol_bc
	rm -f $(SIM)EOSCoeff
	rm -f $(SIM)Comp
	rm -f $(SIM)Comp2
	rm -f $(SIM)suf_comp_bc
	rm -f gmon.out
	rm -f *.d.*   
	rm -f time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
