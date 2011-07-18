input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp adjoint_A.swml adjoint_B.swml; sed -i 's/_A/_B/' adjoint_B.swml; spud-set adjoint_B.swml /timestepping/timestep 0.125
	cp adjoint_A.swml adjoint_C.swml; sed -i 's/_A/_C/' adjoint_C.swml; spud-set adjoint_C.swml /timestepping/timestep 0.0625
	cp adjoint_A.swml adjoint_D.swml; sed -i 's/_A/_D/' adjoint_D.swml; spud-set adjoint_D.swml /timestepping/timestep 0.03125
	cp adjoint_A.swml adjoint_E.swml; sed -i 's/_A/_E/' adjoint_E.swml; spud-set adjoint_E.swml /timestepping/timestep 0.015625

clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc *.pkl 

clean_all:
	rm -f *.vtu *.stat *~*~ *.html *.pyc adjoint_B.swml adjoint_C.swml adjoint_D.swml adjoint_E.swml
	cd src; rm *ele *node *bound *~*~
