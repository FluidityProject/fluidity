PROJECT = rossby

input: clean
	./generate_mesh 25
	./generate_mesh_convergence
	sed s/channel/channel5/ <bdfm1_rossby_wave.swml > rossby5.swml
	sed s/channel/channel10/ <bdfm1_rossby_wave.swml > rossby10.swml
	sed s/channel/channel20/ <bdfm1_rossby_wave.swml > rossby20.swml
	sed s/channel/channel50/ <bdfm1_rossby_wave.swml > rossby50.swml
	sed s/channel/channel100/ <bdfm1_rossby_wave.swml > rossby100.swml
	sed -i s/rossbyx/rossby5/ rossby5.swml
	sed -i s/rossbyx/rossby10/ rossby10.swml
	sed -i s/rossbyx/rossby20/ rossby20.swml
	sed -i s/rossbyx/rossby50/ rossby50.swml
	sed -i s/rossbyx/rossby100/ rossby100.swml

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.edge *.stat *.halo
	rm -f *.msh *.geo
