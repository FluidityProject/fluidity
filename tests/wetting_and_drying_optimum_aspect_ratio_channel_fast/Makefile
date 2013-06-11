input:
	@for f in 1; do cd $$f; ln -sf ../channel.msh .; cd ..; done

clean:
	@rm -f 1/channel{_0.vtu,_1.vtu,.stat} 1/fluidity{.err-0,.log-0} 1/channel.msh 1/matrixdump{,.info}

