run:	input
	mpiexec -n 6 --output-filename log12 ../../bin/shallow_water_mimetic -v 3 W2_12.swml
	mpiexec -n 6 --output-filename log24 ../../bin/shallow_water_mimetic -v 3 W2_24.swml
	mpiexec -n 6 --output-filename log48 ../../bin/shallow_water_mimetic -v 3 W2_48.swml
	mpiexec -n 6 --output-filename log96 ../../bin/shallow_water_mimetic -v 3 W2_96.swml

input: clean
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube12 1 12 --project equiangular
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube24 1 24 --project equiangular
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube48 1 48 --project equiangular
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube96 1 96 --project equiangular
	sed -e 's/{SIMULATION}/W2_12/g' -e 's/{MESH}/cube12/g' -e 's/{TIMESTEP}/7200.0/g' W2.swml > W2_12.swml
	sed -e 's/{SIMULATION}/W2_24/g' -e 's/{MESH}/cube24/g' -e 's/{TIMESTEP}/3600.0/g' W2.swml > W2_24.swml
	sed -e 's/{SIMULATION}/W2_48/g' -e 's/{MESH}/cube48/g' -e 's/{TIMESTEP}/1800.0/g' W2.swml > W2_48.swml
	sed -e 's/{SIMULATION}/W2_96/g' -e 's/{MESH}/cube96/g' -e 's/{TIMESTEP}/900.0/g' W2.swml > W2_96.swml

clean:
	rm -f *.stat *.vtu log *.halo *.node *.ele *.debug *.pyc