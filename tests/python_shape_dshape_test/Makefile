input: clean meshes copy_tests

meshes:
	cd src; make

copy_tests:
	cp test_A.swml test_B.swml; sed -i 's/_A/_B/' test_B.swml
	cp test_A.swml test_C.swml; sed -i 's/_A/_C/' test_C.swml
	cp test_A.swml test_D.swml; sed -i 's/_A/_D/' test_D.swml

clean:
	rm -f *.vtu *.stat src/*.bound src/*.ele src/*.node test_[BCD].swml
