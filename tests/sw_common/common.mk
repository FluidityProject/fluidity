ifndef OUTPUTMESH
OUTPUTMESH=VelocityMesh
endif

sphere_mesh:
	tar -xzvf sw_common/sphere.tgz

fsphere_swml:
	sed -e 's/{MESH}/sphere_ico3/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' sw_common/$(SWML) > fsphere.swml

channel_mesh_convergence:
	sw_common/generate_mesh_convergence

channel_mesh_convergence_quads:
	sw_common/generate_mesh_convergence --quads

rossby_swml:
	sed -e 's/{SIMULATION}/$(SIMULATION)5/g' -e 's/{MESH}/channel5/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' -e 's/{OUTPUTMESH}/$(OUTPUTMESH)/g' sw_common/$(SWML) > $(SIMULATION)5.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)10/g' -e 's/{MESH}/channel10/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' -e 's/{OUTPUTMESH}/$(OUTPUTMESH)/g' sw_common/$(SWML) > $(SIMULATION)10.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)20/g' -e 's/{MESH}/channel20/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' -e 's/{OUTPUTMESH}/$(OUTPUTMESH)/g' sw_common/$(SWML) > $(SIMULATION)20.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)50/g' -e 's/{MESH}/channel50/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' -e 's/{OUTPUTMESH}/$(OUTPUTMESH)/g' sw_common/$(SWML) > $(SIMULATION)50.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)100/g' -e 's/{MESH}/channel100/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' -e 's/{OUTPUTMESH}/$(OUTPUTMESH)/g' sw_common/$(SWML) > $(SIMULATION)100.swml

clean:
	rm -f *.vtu *.stat *.log *.geo *.msh *.ele *.edge *.node *.pyc *.swml
