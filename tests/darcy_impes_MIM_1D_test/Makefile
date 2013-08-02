DIMLMODEL = darcy_transport_immobile_mobile_1D_test

default: input

input: clean
	interval --dx=0.01 --region_ids=region.py 0 0.3 line
	

clean: 
	rm -rf $(DIMLMODEL)_*
	rm -f *txt
	rm -f *.node
	rm -f *.ele
	rm -f *.bound
	rm -f *.vtu
	rm -f *.stat
	rm -f *.log*
	rm -f *.err*
	rm -f matrixdump*
	rm -f src/*.png
