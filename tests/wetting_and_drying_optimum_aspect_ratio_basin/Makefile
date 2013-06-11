input:
	@for f in 1 2 3; do cd $$f; ln -sf ../squarebasin.msh .; cd ..; done

clean:
	@rm -f {1,2,3}/basin{_0.vtu,_1.vtu,.stat} {1,2,3}/fluidity{.err-0,.log-0} {1,2,3}/squarebasin.msh {1,2,3}/matrixdump{,.info}

