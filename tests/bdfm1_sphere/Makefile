default: input

input: clean
	tar -xzvf src/sphere.tgz
	sed s/sphere_ico4/sphere_ico3/ <linear_williamson2_projection.swml > projection3.swml
	sed s/sphere_ico4/sphere_ico4/ <linear_williamson2_projection.swml > projection4.swml
	sed s/sphere_ico4/sphere_ico5/ <linear_williamson2_projection.swml > projection5.swml
	sed s/sphere_ico4/sphere_ico6/ <linear_williamson2_projection.swml > projection6.swml
	sed -i s/LW2/projection3/ projection3.swml
	sed -i s/LW2/projection4/ projection4.swml
	sed -i s/LW2/projection5/ projection5.swml
	sed -i s/LW2/projection6/ projection6.swml
clean:
	rm -f *.vtu *.s *.stat *.log *.msh *.node *.ele *.edge
