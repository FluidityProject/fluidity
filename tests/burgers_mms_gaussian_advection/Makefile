input: clean

copy_tests:
	cp mms_A.bml mms_B.bml; sed -i 's/_A/_B/' mms_B.bml
	cp mms_A.bml mms_C.bml; sed -i 's/_A/_C/' mms_C.bml
	cp mms_A.bml mms_D.bml; sed -i 's/_A/_D/' mms_D.bml
	cp mms_A.bml mms_E.bml; sed -i 's/_A/_E/' mms_E.bml

clean:
	rm -f *.vtu *.stat
