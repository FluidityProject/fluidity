input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp adjoint_A.swml adjoint_B.swml; sed -i 's/_A/_B/' adjoint_B.swml; spud-set adjoint_B.swml /timestepping/timestep 0.25
	cp adjoint_A.swml adjoint_C.swml; sed -i 's/_A/_C/' adjoint_C.swml; spud-set adjoint_C.swml /timestepping/timestep 0.125

clean:
	rm -f *.vtu *.stat
