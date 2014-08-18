DIMLMODEL = Chalcopyrite_dissolution_1D_test

default: input

input: clean
	interval --dx=0.05 --region_ids=region.py 0 0.2 line
	

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
	rm -f src/*.png
