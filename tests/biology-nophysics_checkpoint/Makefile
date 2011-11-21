input: clean
	../../bin/gmsh2triangle src/test1.msh
clean:
	rm -f  *.poly *.s *.d.1 *.detectors
	rm -f test1_checkpoint.stat
	rm -f test1_checkpoint*.vtu

# this cleans everything up so you have to remake the checkpointed data (nodes, etc)
# for the checkpointed run by running the test1.flml
reallyClean:
	rm -f  *.ele *.edge *.face *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors
