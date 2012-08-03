input: clean
	tar -xzvf ../bdfm1_sphere/src/sphere.tgz

	cp galewsky_unperturbed.swml galewsky_unperturbed_3.swml
	cp galewsky_unperturbed.swml galewsky_unperturbed_4.swml
	cp galewsky_unperturbed.swml galewsky_unperturbed_5.swml
	cp galewsky_unperturbed.swml galewsky_unperturbed_6.swml
	spud-set galewsky_unperturbed_3.swml /simulation_name 'sphere_galewsky_unperturbedmesh3'
	spud-set galewsky_unperturbed_4.swml /simulation_name 'sphere_galewsky_unperturbedmesh4'
	spud-set galewsky_unperturbed_5.swml /simulation_name 'sphere_galewsky_unperturbedmesh5'
	spud-set galewsky_unperturbed_6.swml /simulation_name 'sphere_galewsky_unperturbedmesh6'
	sed -i s/sphere_ico4/sphere_ico3/ galewsky_unperturbed_3.swml
	sed -i s/sphere_ico4/sphere_ico4/ galewsky_unperturbed_4.swml
	sed -i s/sphere_ico4/sphere_ico5/ galewsky_unperturbed_5.swml
	sed -i s/sphere_ico4/sphere_ico6/ galewsky_unperturbed_6.swml

	cp galewsky.swml galewsky_3.swml
	cp galewsky.swml galewsky_4.swml
	cp galewsky.swml galewsky_5.swml
	cp galewsky.swml galewsky_6.swml
	spud-set galewsky_3.swml /simulation_name 'sphere_galewskymesh3'
	spud-set galewsky_4.swml /simulation_name 'sphere_galewskymesh4'
	spud-set galewsky_5.swml /simulation_name 'sphere_galewskymesh5'
	spud-set galewsky_6.swml /simulation_name 'sphere_galewskymesh6'
	sed -i s/sphere_ico4/sphere_ico3/ galewsky_3.swml
	sed -i s/sphere_ico4/sphere_ico4/ galewsky_4.swml
	sed -i s/sphere_ico4/sphere_ico5/ galewsky_5.swml
	sed -i s/sphere_ico4/sphere_ico6/ galewsky_6.swml

galewsky_unperturbed_runs: input
	../../bin/shallow_water_mimetic -v 1 galewsky_unperturbed_3.swml &> 3.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_unperturbed_4.swml &> 4.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_unperturbed_5.swml &> 5.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_unperturbed_6.swml &> 6.out &

galewsky_runs: input
	../../bin/shallow_water_mimetic -v 1 galewsky_3.swml &> 3.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_4.swml &> 4.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_5.swml &> 5.out &
	../../bin/shallow_water_mimetic -v 1 galewsky_6.swml &> 6.out &

clean:
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f channel*.ele channel*.face chanel*.node channel*.edge *.stat *.halo
	rm -f *.msh *.geo
