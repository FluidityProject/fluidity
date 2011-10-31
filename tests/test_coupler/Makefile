preprocess:
	@echo ********Generating the mesh using gmsh in 2d:
	gmsh -2 ./domain.geo
	@echo ********Converting the gmsh mesh to triangle format:
	../../bin/gmsh2triangle --2d ./domain.msh

	@echo ********Generating the mesh for OASIS coupler
	gmsh -2 ./domain_struct.geo
	@echo ********Converting the gmsh mesh to triangle format:
	../../bin/gmsh2triangle --2d ./domain_struct.msh

clean: 
	rm -f  *.msh *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump*

run:
	@echo **********Calling test_coupler w/out logging:
	../../bin/test_coupler test_coupler.cml

input: clean preprocess