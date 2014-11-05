#include "petsc.h"
#include "petscsf.h"

PetscErrorCode dmplex_get_mesh_connectivity(DM plex, PetscInt nnodes, PetscInt loc, IS *rnbrCells, IS *rnbrVertices, PetscInt *ndglno);
PetscErrorCode dmplex_get_num_surface_facets(DM plex, const char labelname[], PetscInt *nfacets);
PetscErrorCode dmplex_get_face_connectivity(DM plex, const char labelname[], PetscInt nfacets, PetscInt sloc, PetscInt *sndglno, PetscInt *boundary_ids);
