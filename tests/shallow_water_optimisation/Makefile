input: meshes copy_tests clean_input

meshes:
	cd src; make

copy_tests:
	cp adjoint_template.swml adjoint_A.swml; sed -i 's/_template/_A/' adjoint_A.swml
	cp adjoint_template.swml adjoint_B.swml; sed -i 's/_template/_B/' adjoint_B.swml
	cp adjoint_template.swml adjoint_C.swml; sed -i 's/_template/_C/' adjoint_C.swml
	cp adjoint_template.swml adjoint_D.swml; sed -i 's/_template/_D/' adjoint_D.swml
	cp adjoint_template.swml adjoint_E.swml; sed -i 's/_template/_E/' adjoint_E.swml
	spud-set adjoint_A.swml /timestepping/timestep 0.25
	spud-set adjoint_B.swml /timestepping/timestep 0.125
	spud-set adjoint_C.swml /timestepping/timestep 0.0625
	spud-set adjoint_D.swml /timestepping/timestep 0.03125
	spud-set adjoint_E.swml /timestepping/timestep 0.015625
	cp optimiser_template.oml optimiser_A.oml; sed -i 's/_template/_A/' optimiser_A.oml
	cp optimiser_template.oml optimiser_B.oml; sed -i 's/_template/_B/' optimiser_B.oml
	cp optimiser_template.oml optimiser_C.oml; sed -i 's/_template/_C/' optimiser_C.oml
	cp optimiser_template.oml optimiser_D.oml; sed -i 's/_template/_D/' optimiser_D.oml
	cp optimiser_template.oml optimiser_E.oml; sed -i 's/_template/_E/' optimiser_E.oml

clean_input:
	rm -f *.vtu *.stat *~*~ *.html *.pyc *.pkl

clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc adjoint_?.swml optimiser_?.oml *.pkl
	rm -Rf opt_iter_*
	cd src; rm -f *ele *node *bound *~*~
