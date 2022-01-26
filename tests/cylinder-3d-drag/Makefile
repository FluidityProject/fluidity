input: clean
	cp src/cylinder.msh .

clean:
	rm -f *.msh *.poly *.s *.d.* *.stat *.s *.log-0 *.err-0
	find . -name '*.vtu' | xargs -n 20 rm -f
