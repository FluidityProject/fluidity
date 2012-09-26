input: meshes copy_tests clean_input

meshes:
	cd src; make

copy_tests:
	cp adjoint_template.swml adjoint_A.swml
	cp adjoint_template.swml adjoint_B.swml
	cp adjoint_template.swml adjoint_C.swml
	sed -i 's/_template/_A/' adjoint_A.swml
	sed -i 's/_template/_B/' adjoint_B.swml
	sed -i 's/_template/_C/' adjoint_C.swml
	spud-set adjoint_A.swml /timestepping/timestep 0.25
	spud-set adjoint_B.swml /timestepping/timestep 0.125
	spud-set adjoint_C.swml /timestepping/timestep 0.0625

clean_input:
	rm -f *.vtu *.stat *~*~ *.html *.pyc *.pkl

clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc adjoint_?.swml *.pkl
	cd src; rm -f *ele *node *bound *~*~ *.edge *.msh
