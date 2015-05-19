DIMLMODEL = leaching_2phase_heat_transfer_1D_test

default: input

input: clean
	 interval --dx=2.0 --region_ids=./region.py 0 400.0 line
	

clean: 
	rm -rf $(DIMLMODEL)_*
	rm -f *.node
	rm -f *.ele
	rm -f *.bound
	rm -f *.vtu
	rm -f *.stat
	rm -f *.log*
	rm -f *.err*
	rm -f matrixdump*

