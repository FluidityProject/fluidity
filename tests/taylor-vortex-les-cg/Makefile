
input: clean
	gmsh -3 -optimize -o MMS_$(N).msh Stokes-test-files/MMS_$(N).geo
	../../bin/gmsh2triangle -i MMS_$(N).msh
	../../bin/periodise taylor-green-vortex-3d-periodic-$(N).flml
	mpiexec -n $(P) ../../bin/flredecomp -i 1 -o $(P) taylor-green-vortex-3d-periodic-$(N)_periodised taylor-green-vortex-3d-periodic-$(N)-$(P)p
	#../../bin/create_aligned_mesh cube_0_2pi_12 6.283185307179586 6.283185307179586 6.283185307179586 12 12 12
	cd ..

clean:
	#rm -f *.ele *.face *.node
