preprocess:

run:
	@echo **********Running simulation
	mpiexec -n 16 flredecomp -i 1 -o 16 -v -l flow_past_a_square flow_past_a_square_flredecomp
	mpiexec -n 16 fluidity flow_past_a_square_flredecomp.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *halo
	rm -f *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf flow_past_a_square_*
	rm -f nohup.out
	rm -rf *flredecomp*
