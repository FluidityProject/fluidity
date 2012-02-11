
default: run

run: input
	../../bin/fluidity -v 3 -l missingdata.flml 

%.flml: valid.flml
	@sed -e 's/valid\.nc/$*.nc/' $^ > $@

input: clean
	./createnetcdf.py
	$(MAKE) missingdata.flml missingdimension.flml missingvariable.flml incorrectdimension.flml incorrectvariable.flml
	../../bin/create_aligned_mesh mesh 2.0 2.0 1.0 10 10 1 -- -1.0 -1.0 -1.0

clean:
	rm -f *.{nc,pyc,ele,node,face,vtu,stat,s,err} matrixdump* fluidity.*
	rm -f missing{data,dimension,variable}.flml incorrect{dimension,variable}.flml

