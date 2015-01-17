input: MMS_A.msh MMS_B.msh MMS_C.msh MMS_D.msh MMS_E.msh


MMS_%.msh: src/MMS_%.geo
	gmsh -2 $< -o $@

clean:
	rm -rf output* fine* diff* *.msh fluidity.* reference*
