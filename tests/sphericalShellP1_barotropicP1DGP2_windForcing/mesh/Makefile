# Definition of directory locations (Overriden if recursive make is used from the parent
# directory, with values defined therein):
MESH_FILE_PATTERN=sphericalShell
FLUIDITY=../../../
BIN=$(FLUIDITY)bin
SCHEMAS=$(FLUIDITY)schemas

GMSH_GEO_FILE=$(MESH_FILE_PATTERN).geo
GMSH_MSH_FILE=$(MESH_FILE_PATTERN).msh

all: mesh

mesh:$(GMSH_MSH_FILE)
	@echo "**********Converting the gmsh mesh (ASCII) to triangle format:"
	$(BIN)/gmsh2triangle -s $<

$(GMSH_MSH_FILE):$(GMSH_GEO_FILE)
	gmsh -2 $<

clean:
	rm -f *.edge *.ele *.node *.halo $(GMSH_MSH_FILE)
	rm -f *[0-9].egde *[0-9].ele *[0-9].node *[0-9].halo
