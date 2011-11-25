input: clean copy_tests

copy_tests:
	cp mms_A.bml mms_B.bml; sed -i 's/_A/_B/' mms_B.bml
	cp mms_A.bml mms_C.bml; sed -i 's/_A/_C/' mms_C.bml
	cp mms_A.bml mms_D.bml; sed -i 's/_A/_D/' mms_D.bml
	cp mms_A.bml mms_E.bml; sed -i 's/_A/_E/' mms_E.bml
	spud-set mms_A.bml /timestepping/timestep 2.0
	spud-set mms_B.bml /timestepping/timestep 1.0
	spud-set mms_C.bml /timestepping/timestep 0.5
	spud-set mms_D.bml /timestepping/timestep 0.25
	spud-set mms_E.bml /timestepping/timestep 0.125

clean:
	rm -f *.vtu *.stat
