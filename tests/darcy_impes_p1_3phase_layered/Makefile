include ../darcy_impes_common/Depend.mk

SIM = darcy_impes_p1_3phase_layered

input: clean
	cp mesh_data/line*.msh .
	cp mesh_data/square*.msh .
	cp mesh_data/cube*.msh .

clean:
	rm -f *.ele
	rm -f *.node
	rm -f *.bound
	rm -f *.edge
	rm -f *.face
	rm -f *.vtu 
	rm -f *.pvtu 
	rm -f *.s 
	rm -f *.stat 
	rm -f *.log-0 
	rm -f *.err-0 
	rm -f *.msh
	rm -f *.halo
	rm -f $(SIM)*.err*
	rm -f $(SIM)*.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f *.convergence 
	rm -f error_*.txt
