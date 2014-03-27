input: clean
	tar -xzvf ../sw_common/sphere.tgz
	cp sphere_williamson2_new.swml sphere_williamson2_new3.swml
	cp sphere_williamson2_new.swml sphere_williamson2_new4.swml
	cp sphere_williamson2_new.swml sphere_williamson2_new5.swml
	cp sphere_williamson2_new.swml sphere_williamson2_new6.swml
	spud-set sphere_williamson2_new3.swml /simulation_name 'sphere_williamson2mesh3'
	spud-set sphere_williamson2_new4.swml /simulation_name 'sphere_williamson2mesh4'
	spud-set sphere_williamson2_new5.swml /simulation_name 'sphere_williamson2mesh5'
	spud-set sphere_williamson2_new6.swml /simulation_name 'sphere_williamson2mesh6'
	sed -i s/sphere_ico4/sphere_ico3/ sphere_williamson2_new3.swml
	sed -i s/sphere_ico4/sphere_ico4/ sphere_williamson2_new4.swml
	sed -i s/sphere_ico4/sphere_ico5/ sphere_williamson2_new5.swml
	sed -i s/sphere_ico4/sphere_ico6/ sphere_williamson2_new6.swml

	spud-set sphere_williamson2_new5.swml /timestepping/timestep 900
	spud-set sphere_williamson2_new6.swml /timestepping/timestep 450
clean:
	rm -f  *.vtu *.stat
	rm -f *.ele *.node
	rm -f *.out
	rm -f sphere_williamson2_new3.swml
	rm -f sphere_williamson2_new4.swml
	rm -f sphere_williamson2_new5.swml
	rm -f sphere_williamson2_new6.swml

sw2_convergence_runs: input
	../../bin/shallow_water_mimetic -v 1 sphere_williamson2_new3.swml &> 3.out &
	../../bin/shallow_water_mimetic -v 1 sphere_williamson2_new4.swml &> 4.out &
	../../bin/shallow_water_mimetic -v 1 sphere_williamson2_new5.swml &> 5.out &
	../../bin/shallow_water_mimetic -v 1 sphere_williamson2_new6.swml &> 6.out &