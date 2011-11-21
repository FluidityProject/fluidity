PROJECT = channel_parallel

input: clean
	gmsh -3 -o channel.msh src/channel.geo
	../../bin/gmsh2triangle channel.msh
	../../bin/periodise channel_parallel.flml
	mpiexec -np 4 ../../bin/flredecomp -i 1 -o 4 channel_parallel_periodised channel_parallel_periodised_flredecomp

clean:
	rm -f $(PROJECT).stat
	rm -f *vtu
	rm -f fluidity*
	rm -f *.ele *.face *.node 
	rm -f *.msh
	rm -f *_periodised.flml
	rm -f *_flredecomp.flml
	rm -f *.halo
	rm -rf channel_parallel_[0-4]*
