input: clean copy_tests

copy_tests:
	cp mms_template.bml mms_A.bml; sed -i 's/_template/_A/' mms_A.bml
	cp mms_template.bml mms_B.bml; sed -i 's/_template/_B/' mms_B.bml
	cp mms_template.bml mms_C.bml; sed -i 's/_template/_C/' mms_C.bml
	cp mms_template.bml mms_D.bml; sed -i 's/_template/_D/' mms_D.bml
	cp mms_template.bml mms_E.bml; sed -i 's/_template/_E/' mms_E.bml
	spud-set mms_A.bml /timestepping/timestep 2.0
	spud-set mms_B.bml /timestepping/timestep 1.0
	spud-set mms_C.bml /timestepping/timestep 0.5
	spud-set mms_D.bml /timestepping/timestep 0.25
	spud-set mms_E.bml /timestepping/timestep 0.125

clean:
	rm -f *.vtu *.stat mms_?.bml
