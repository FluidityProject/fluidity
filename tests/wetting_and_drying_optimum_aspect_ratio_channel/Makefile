input:
	@for f in 1 2 3 4 5; do cd $$f; ln -sf ../channel.msh .; cd ..; done

clean:
	@rm -f {1,2,3,4,5}/channel{_0.vtu,_1.vtu,.stat} {1,2,3,4,5}/fluidity{.err-0,.log-0} {1,2,3,4,5}/channel.msh {1,2,3,4,5}/matrixdump{,.info}

