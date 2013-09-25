default: input

input: clean
	tar -xzvf ../sw_common/sphere.tgz
	sed s/spheremesh/sphere_ico3/ <sw_fsphere.swml > fsphere3.swml
	sed s/spheremesh/sphere_ico4/ <sw_fsphere.swml > fsphere4.swml
	sed s/spheremesh/sphere_ico5/ <sw_fsphere.swml > fsphere5.swml
	sed s/spheremesh/sphere_ico6/ <sw_fsphere.swml > fsphere6.swml
	sed -i s/fsphere/fsphere3/ fsphere3.swml
	sed -i s/fsphere/fsphere4/ fsphere4.swml
	sed -i s/fsphere/fsphere5/ fsphere5.swml
	sed -i s/fsphere/fsphere6/ fsphere6.swml
clean:
	rm -f *.vtu *.stat *.log sphere*.ele sphere*.node fsphere*.swml
