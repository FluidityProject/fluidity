
default: run

run: input
	../../bin/fluidity -v 3 -l netcdf_read_free_surface.flml 

missing%.flml: valid.flml
	sed -e 's/valid\.nc/missing$*.nc/' $^ > $@

input: clean
	./createnetcdf.py
	$(MAKE) missingdata.flml missingdimension.flml missingvariable.flml
	../../bin/create_aligned_mesh mesh 2.0 2.0 1.0 10 10 1 -- -1.0 -1.0 -1.0

clean:
	rm -f *.nc *.pyc *.ele *.node *.face *.vtu *.stat *.s matrixdump* fluidity.*
	rm -f missing{data,dimension,variable}.flml 

