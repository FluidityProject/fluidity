input: clean
	#gmsh -3 -o cube-12.msh src/cube_12.geo
	#gmsh -3 -o cube-24.msh src/cube_24.geo
	#gmsh -3 -o cube-48.msh src/cube_48.geo
	#../../bin/gmsh2triangle cube-12.msh
	#../../bin/gmsh2triangle cube-24.msh
	#../../bin/gmsh2triangle cube-48.msh
	../../bin/create_aligned_mesh cube_0_2pi_12 6.283185307179586 6.283185307179586 6.283185307179586 12 12 12
	../../bin/create_aligned_mesh cube_0_2pi_24 6.283185307179586 6.283185307179586 6.283185307179586 24 24 24
	../../bin/create_aligned_mesh cube_0_2pi_48 6.283185307179586 6.283185307179586 6.283185307179586 48 48 48
	#../../bin/fldecomp -n 2 -f cube_0_2pi_24
	#../../bin/fldecomp -n 4 -f cube_0_2pi_48

clean:
	#rm -f *.ele *.face *.node
