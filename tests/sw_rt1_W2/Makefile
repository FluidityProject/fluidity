input: clean
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube12 1 12 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube12 6
	mv cube12_serial.node cube12.node
	mv cube12_serial.ele cube12.ele
	rm cube12_*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube24 1 24 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube24 6
	mv cube24_serial.node cube24.node
	mv cube24_serial.ele cube24.ele
	rm cube24_*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube48 1 48 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube48 6
	mv cube48_serial.node cube48.node
	mv cube48_serial.ele cube48.ele
	rm cube48_*
	mpiexec -n 6 ../../tools/generate_cubed_sphere_mesh.py cube96 1 96 --project equiangular
	../../tools/make_serial_cubed_sphere_mesh.py cube96 6
	mv cube96_serial.node cube96.node
	mv cube96_serial.ele cube96.ele
	rm cube96_*
	rm *.debug
	cp sphere_rt1_W2.swml sphere_rt1_W2_12.swml 
	cp sphere_rt1_W2.swml sphere_rt1_W2_24.swml 
	cp sphere_rt1_W2.swml sphere_rt1_W2_48.swml 
	cp sphere_rt1_W2.swml sphere_rt1_W2_96.swml 
	spud-set sphere_rt1_W2.swml /simulation_name 'rt1_W2_12'
	spud-set sphere_rt1_W2.swml /simulation_name 'rt1_W2_24'
	spud-set sphere_rt1_W2.swml /simulation_name 'rt1_W2_48'
	spud-set sphere_rt1_W2.swml /simulation_name 'rt1_W2_96'
	sed -i s/cube_serial/cube12/ sphere_rt1_W2_12.swml
	sed -i s/cube_serial/cube24/ sphere_rt1_W2_24.swml
	sed -i s/cube_serial/cube48/ sphere_rt1_W2_48.swml
	sed -i s/cube_serial/cube96/ sphere_rt1_W2_96.swml

	spud-set sphere_rt1_W2_48.swml /timestepping/timestep 900
	spud-set sphere_rt1_W2_96.swml /timestepping/timestep 450


clean:
	rm -f  *.vtu *.stat
	rm -f *.ele *.node *.halo *.debug
	rm -f *.out
	rm -f sphere_rt1_W2_12.swml
	rm -f sphere_rt1_W2_24.swml
	rm -f sphere_rt1_W2_48.swml
	rm -f sphere_rt1_W2_96.swml

sw2_convergence_runs: input
	../../bin/shallow_water_mimetic -v 3 sphere_rt1_W2_12.swml &> 12.out &
	../../bin/shallow_water_mimetic -v 3 sphere_rt1_W2_24.swml &> 24.out &
	../../bin/shallow_water_mimetic -v 3 sphere_rt1_W2_48.swml &> 48.out &
	../../bin/shallow_water_mimetic -v 3 sphere_rt1_W2_96.swml &> 96.out &
