input: clean
	tar -xzvf ../bdfm1_sphere/src/sphere.tgz

	cp williamson5.swml williamson5_3.swml
	cp williamson5.swml williamson5_4.swml
	cp williamson5.swml williamson5_5.swml
	cp williamson5.swml williamson5_6.swml
	spud-set williamson5_3.swml /simulation_name 'sphere_williamson5mesh3'
	spud-set williamson5_4.swml /simulation_name 'sphere_williamson5mesh4'
	spud-set williamson5_5.swml /simulation_name 'sphere_williamson5mesh5'
	spud-set williamson5_6.swml /simulation_name 'sphere_williamson5mesh6'
	sed -i s/sphere_ico4/sphere_ico3/ williamson5_3.swml
	sed -i s/sphere_ico4/sphere_ico4/ williamson5_4.swml
	sed -i s/sphere_ico4/sphere_ico5/ williamson5_5.swml
	sed -i s/sphere_ico4/sphere_ico6/ williamson5_6.swml

	spud-set williamson5_6.swml /timestepping/timestep 900.0

sw2_convergence_runs: input
	../../bin/shallow_water_mimetic -v 1 williamson5_3.swml &> 3.out &
	../../bin/shallow_water_mimetic -v 1 williamson5_4.swml &> 4.out &
	../../bin/shallow_water_mimetic -v 1 williamson5_5.swml &> 5.out &
	../../bin/shallow_water_mimetic -v 1 williamson5_6.swml &> 6.out &

sw2_6: input
	../../bin/shallow_water_mimetic -v 1 williamson5_6.swml &> 6.out &

clean:
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f channel*.ele channel*.face chanel*.node channel*.edge *.stat *.halo
	rm -f *.msh *.geo
