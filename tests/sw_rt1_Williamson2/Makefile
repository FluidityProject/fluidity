input: clean
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube12 1 12 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube12 6
	rm cube12_0.* cube12_1.* cube12_2.* cube12_3.* cube12_4.* cube12_5.*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube24 1 24 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube24 6
	rm cube24_0.* cube24_1.* cube24_2.* cube24_3.* cube24_4.* cube24_5.*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube48 1 48 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube48 6
	rm cube48_0.* cube48_1.* cube48_2.* cube48_3.* cube48_4.* cube48_5.*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube96 1 96 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube96 6
	rm cube96_0.* cube96_1.* cube96_2.* cube96_3.* cube96_4.* cube96_5.*
	sed -e 's/{SIMULATION}/W2_12/g' -e 's/{MESH}/cube12_serial/g' -e 's/{TIMESTEP}/7200.0/g' W2.swml > W2_12.swml
	sed -e 's/{SIMULATION}/W2_24/g' -e 's/{MESH}/cube24_serial/g' -e 's/{TIMESTEP}/3600.0/g' W2.swml > W2_24.swml
	sed -e 's/{SIMULATION}/W2_48/g' -e 's/{MESH}/cube48_serial/g' -e 's/{TIMESTEP}/1800.0/g' W2.swml > W2_48.swml
	sed -e 's/{SIMULATION}/W2_96/g' -e 's/{MESH}/cube96_serial/g' -e 's/{TIMESTEP}/900.0/g' W2.swml > W2_96.swml

clean:
	rm -f *.stat *.vtu log *.halo *.node *.ele *.debug *.pyc