SHELL = sh

SIM = rad_1d_1grp_hom_vacbc_eig

# options for 1d mesh interval script for different mesh resolutions
LEFT_X    = 0.0
RIGHT_X   = 4.513502

MESH_SIZE_A = 2.256751
MESH_SIZE_B = 1.1283755
MESH_SIZE_C = 0.56418775
MESH_SIZE_D = 0.282093875
MESH_SIZE_E = 0.141046937
MESH_SIZE_F = 0.070523469
MESH_SIZE_G = 0.035261734
MESH_SIZE_H = 0.017630867
MESH_SIZE_I = 0.008815434
MESH_SIZE_J = 0.004407717
MESH_SIZE_K = 0.002203858

input: clean
	@echo "Running INTERVAL script for $(SIM) for each mesh resolution"
	../../scripts/interval --dx=$(MESH_SIZE_A) $(LEFT_X) $(RIGHT_X) $(SIM)_A
	../../scripts/interval --dx=$(MESH_SIZE_B) $(LEFT_X) $(RIGHT_X) $(SIM)_B
	../../scripts/interval --dx=$(MESH_SIZE_C) $(LEFT_X) $(RIGHT_X) $(SIM)_C
	../../scripts/interval --dx=$(MESH_SIZE_D) $(LEFT_X) $(RIGHT_X) $(SIM)_D
	../../scripts/interval --dx=$(MESH_SIZE_E) $(LEFT_X) $(RIGHT_X) $(SIM)_E
	../../scripts/interval --dx=$(MESH_SIZE_F) $(LEFT_X) $(RIGHT_X) $(SIM)_F
	../../scripts/interval --dx=$(MESH_SIZE_G) $(LEFT_X) $(RIGHT_X) $(SIM)_G
	../../scripts/interval --dx=$(MESH_SIZE_H) $(LEFT_X) $(RIGHT_X) $(SIM)_H
	../../scripts/interval --dx=$(MESH_SIZE_I) $(LEFT_X) $(RIGHT_X) $(SIM)_I
	../../scripts/interval --dx=$(MESH_SIZE_J) $(LEFT_X) $(RIGHT_X) $(SIM)_J
	../../scripts/interval --dx=$(MESH_SIZE_K) $(LEFT_X) $(RIGHT_X) $(SIM)_K

clean:
	rm -f $(SIM)*.ele
	rm -f $(SIM)*.node
	rm -f $(SIM)*.bound
	rm -f $(SIM)*.edge
	rm -f $(SIM)*.face
	rm -f $(SIM)*.vtu 
	rm -f $(SIM)*.pvtu 
	rm -f $(SIM)*.s 
	rm -f $(SIM)*.d.* 
	rm -f $(SIM)*.stat 
	rm -f $(SIM)*.log-0 
	rm -f $(SIM)*.err-0 
	rm -f $(SIM)*.msh
	rm -f $(SIM)*.halo
	rm -f fluidity.err*
	rm -f fluidity.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
