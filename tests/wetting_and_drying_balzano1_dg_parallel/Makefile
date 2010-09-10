all: input
	fluidity  -v5 balzano1_2plus1.flml &> output.txt

input: clean
	cd mesh; ./makemesh_parallel.sh

clean:
	rm -f  *.ele *.edge *.node *.poly *vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.png *.dat
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.halo 
config:
	diamond -s ../../tools/fluidity_options.rng  balzano1_2plus1.flml 

