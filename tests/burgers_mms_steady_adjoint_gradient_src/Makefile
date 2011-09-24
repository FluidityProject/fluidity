input: clean copy_tests

copy_tests:
	cp mms_A.bml mms_B.bml; sed -i 's/_A/_B/' mms_B.bml
	cp mms_A.bml mms_C.bml; sed -i 's/_A/_C/' mms_C.bml
	cp mms_A.bml mms_D.bml; sed -i 's/_A/_D/' mms_D.bml
	cp mms_A.bml mms_E.bml; sed -i 's/_A/_E/' mms_E.bml
	cp op_A.oml op_B.oml; sed -i 's/_A/_B/' op_B.oml
	cp op_A.oml op_C.oml; sed -i 's/_A/_C/' op_C.oml
	cp op_A.oml op_D.oml; sed -i 's/_A/_D/' op_D.oml
	cp op_A.oml op_E.oml; sed -i 's/_A/_E/' op_E.oml

clean:
	rm -f *.vtu *.stat *.pkl *tmp.[ob]ml burgers.py
