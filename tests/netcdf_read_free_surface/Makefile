
default: run

run: input
	../../bin/fluidity -v 3 -l netcdf_read_free_surface.flml 

input: clean
	./createnetcdf.py
	../../bin/create_aligned_mesh mesh 2.0 2.0 1.0 10 10 1 -- -1.0 -1.0 -1.0

clean:
	rm -f height.nc *.pyc *.ele *.node *.face *.vtu *.stat *.s matrixdump* fluidity.*

