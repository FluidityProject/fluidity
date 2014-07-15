input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp adjoint_template.swml adjoint_A.swml; sed -i 's/_template/_A/' adjoint_A.swml; spud-set adjoint_A.swml /timestepping/timestep 0.125
	cp adjoint_template.swml adjoint_B.swml; sed -i 's/_template/_B/' adjoint_B.swml; spud-set adjoint_B.swml /timestepping/timestep 0.0625
	cp adjoint_template.swml adjoint_C.swml; sed -i 's/_template/_C/' adjoint_C.swml; spud-set adjoint_C.swml /timestepping/timestep 0.03125

clean:
	rm -f *.vtu *.stat src/*.bound src/*.ele src/*.node adjoint_?.swml
