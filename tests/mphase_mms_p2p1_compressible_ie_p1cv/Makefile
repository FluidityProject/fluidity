input:	clean
	./copy_script

run:
	@echo **********Running simulations
	fluidity MMS_A.flml
	fluidity MMS_B.flml
	fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* matrix* *~ mphase_mms_p2p1_compressible_ie_p1cv.py *.pyc fluidity.*
