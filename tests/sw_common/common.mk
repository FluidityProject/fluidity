mesh:
	tar -xzvf sw_common/sphere.tgz

swml:
	sed -e 's/{SIMULATION}/$(SIMULATION)_3/g' -e 's/{MESH}/sphere_ico3/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' sw_common/$(SWML) > $(SIMULATION)_3.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)_4/g' -e 's/{MESH}/sphere_ico4/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' sw_common/$(SWML) > $(SIMULATION)_4.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)_5/g' -e 's/{MESH}/sphere_ico5/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' sw_common/$(SWML) > $(SIMULATION)_5.swml
	sed -e 's/{SIMULATION}/$(SIMULATION)_6/g' -e 's/{MESH}/sphere_ico6/g' -e 's/{PDEGREE}/$(PDEGREE)/g' -e 's/{UDEGREE}/$(UDEGREE)/g' -e 's/{VDEGREE}/$(VDEGREE)/g' -e 's/{CONSTRAINT}/$(CONSTRAINT)/g' sw_common/$(SWML) > $(SIMULATION)_6.swml

clean:
	rm -f *.vtu *.stat *.log sphere*.ele sphere*.node *.pyc *.swml
