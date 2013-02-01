input: clean copy_tests

copy_tests:
	cp mms_template.bml mms_A.bml; sed -i 's/_template/_A/' mms_A.bml
	cp mms_template.bml mms_B.bml; sed -i 's/_template/_B/' mms_B.bml
	cp mms_template.bml mms_C.bml; sed -i 's/_template/_C/' mms_C.bml
	cp mms_template.bml mms_D.bml; sed -i 's/_template/_D/' mms_D.bml
	cp mms_template.bml mms_E.bml; sed -i 's/_template/_E/' mms_E.bml
	cp op_template.oml op_A.oml; sed -i 's/_template/_A/' op_A.oml
	cp op_template.oml op_B.oml; sed -i 's/_template/_B/' op_B.oml
	cp op_template.oml op_C.oml; sed -i 's/_template/_C/' op_C.oml
	cp op_template.oml op_D.oml; sed -i 's/_template/_D/' op_D.oml
	cp op_template.oml op_E.oml; sed -i 's/_template/_E/' op_E.oml

clean:
	rm -f *.vtu *.stat *.pkl *tmp.oml mms_?.bml op_?.oml burgers.py
