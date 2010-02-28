default: input

input: clean
	tar -xvjf src/square.tar.bz2
	tar -xvjf src/target_meshes.tar.bz2

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat *.log src/*.msh
	rm -f geostrophic_interpolation_new*vtu geostrophic_interpolation_old*vtu
