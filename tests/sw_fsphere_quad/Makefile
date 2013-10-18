input: clean
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube 1 10 --debug --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube 6

clean:
	rm -f *.stat *.vtu log *.halo *.node *.ele *.debug *.pyc