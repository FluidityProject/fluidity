ROOT = /d/src/fluidity_mixingdg
# ${ROOT}/scripts/gmsh_mesh_transform True '(x, 0.001 * y, z)' square.msh

input: clean
	gmsh -2 src/square.geo -o square.msh
	${ROOT}/scripts/gmsh2triangle --2d square.msh

clean:
	rm -f  *.ele *.edge *.node *.msh *.vtu *.stat matrixdump* *.node.bak fluidity.*

run: input
	${ROOT}/bin/fluidity -v3 -l mix.flml

diamond:
	diamond -s ${ROOT}/tools/fluidity_options.rng mix.flml

#diff .{./mix$(expr $(basename $PWD | sed -e "s/^[^0-9]*\([0-9]*\).*$/\1/") - 1)*,}/mix.flml

view:
	mayavi -z bandsurfMesh.mv

view2:
	mayavi -d mix_0.vtu -m BandedSurfaceMap

viewquick:
	vtu2img.sh -la -ng && img2txt png/latest.png

fluidity:
	cd ${ROOT} && spud-preprocess tools/fluidity_options.rnc && make

all: fluidity run viewquick

