input: clean
	gmsh -3 src/channel.geo
	cp src/channel.msh .

clean:
	rm -rf *.stat *.vtu channel.node channel.ele channel.face *.convergence *checkpoint* fluidity.log-0 fluidity.err-0 src/channel.msh \
	matrixdump matrixdump.info
