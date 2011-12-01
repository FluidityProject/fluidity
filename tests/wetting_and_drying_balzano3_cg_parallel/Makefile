all: input
	fluidity  -v5 balzano1_2plus1.flml &> output.txt

input: clean
	cd mesh; ./makemesh_parallel.sh

clean:
	rm -f  *.ele *.edge *.node *.poly *vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.png *.dat
	rm -rf balzano1_2plus1_? balzano1_2plus1_?? balzano1_2plus1_???
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors matrixdump* *.halo 
config:
	diamond -s ../../tools/fluidity_options.rng  balzano1_2plus1.flml 

