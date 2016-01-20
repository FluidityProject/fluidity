default: input

input: clean
	tar -xvjf src/basin.tar.bz2
	../../bin/triangle2gmsh basin
	tar -xvjf src/first_timestep_adapted_mesh.tar.bz2
	../../bin/triangle2gmsh first_timestep_adapted_mesh

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat src/*.msh
	rm -f gmsh.log gmsh_err.log
clean-run:
	rm -f geostrophic_interpolation_new*vtu geostrophic_interpolation_old*vtu p_decomp_new*vtu p_decomp_old*vtu
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-*
