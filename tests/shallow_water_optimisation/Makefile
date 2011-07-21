input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp adjoint_A.swml adjoint_B.swml; sed -i 's/_A/_B/' adjoint_B.swml; spud-set adjoint_B.swml /timestepping/timestep 0.125
	cp adjoint_A.swml adjoint_C.swml; sed -i 's/_A/_C/' adjoint_C.swml; spud-set adjoint_C.swml /timestepping/timestep 0.0625
	cp adjoint_A.swml adjoint_D.swml; sed -i 's/_A/_D/' adjoint_D.swml; spud-set adjoint_D.swml /timestepping/timestep 0.03125
	cp adjoint_A.swml adjoint_E.swml; sed -i 's/_A/_E/' adjoint_E.swml; spud-set adjoint_E.swml /timestepping/timestep 0.015625
	cp optimiser_A.oml optimiser_B.oml; sed -i 's/_A/_B/' optimiser_B.oml
	cp optimiser_A.oml optimiser_C.oml; sed -i 's/_A/_C/' optimiser_C.oml
	cp optimiser_A.oml optimiser_D.oml; sed -i 's/_A/_D/' optimiser_D.oml
	cp optimiser_A.oml optimiser_E.oml; sed -i 's/_A/_E/' optimiser_E.oml

clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc *.pkl 

clean_all:
	rm -f *.vtu *.stat *~*~ *.html *.pyc adjoint_B.swml adjoint_C.swml adjoint_D.swml adjoint_E.swml optimiser_B.oml optimiser_C.oml optimiser_D.oml optimiser_E.oml 
	rm -R opt_iter_*
	cd src; rm *ele *node *bound *~*~
