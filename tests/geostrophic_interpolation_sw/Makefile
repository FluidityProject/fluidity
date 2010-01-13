default: input

input: clean
	gmsh -2 src/basin.geo > gmsh.log 2> gmsh_err.log
	../../scripts/gmsh2triangle -2 src/basin.msh
	tar -xvjf src/first_timestep_adapted_mesh.tar.bz2

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat *.log *.msh
	rm -f geostrophic_interpolation_new*vtu geostrophic_interpolation_old*vtu
