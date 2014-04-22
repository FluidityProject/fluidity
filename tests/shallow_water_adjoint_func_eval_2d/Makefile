input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp adjoint_template.swml adjoint_A.swml; sed -i 's/_template/_A/' adjoint_A.swml; spud-set adjoint_A.swml /timestepping/timestep 0.5
	cp adjoint_template.swml adjoint_B.swml; sed -i 's/_template/_B/' adjoint_B.swml; spud-set adjoint_B.swml /timestepping/timestep 0.25
	cp adjoint_template.swml adjoint_C.swml; sed -i 's/_template/_C/' adjoint_C.swml; spud-set adjoint_C.swml /timestepping/timestep 0.125

clean:
	rm -f *.vtu *.stat adjoint_?.swml
