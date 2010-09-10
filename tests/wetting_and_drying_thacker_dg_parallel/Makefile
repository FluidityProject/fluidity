all: input
	mpiexec -n 4 fluidity  -v5 thacker.flml &> output.txt

input: clean
	cd mesh; ./makemesh_parallel.sh

clean:
	rm -f  *.ele *.edge *.node *.poly *vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.dat 
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.halo 
config:
	diamond -s ../../tools/fluidity_options.rng  thacker.flml 

