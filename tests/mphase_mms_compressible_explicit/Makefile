input:	clean
	./copy_script

run:
	@echo **********Running simulations
	../../bin/fluidity MMS_A.flml
	../../bin/fluidity MMS_B.flml

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* matrix* *~ mphase_mms_compressible_explicit.py *.pyc fluidity.*
