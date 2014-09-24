#include "petsc.h"

PetscErrorCode dmplex_get_mesh_connectivity(DM plex, PetscInt nnodes, PetscInt loc, PetscInt *ndglno);
PetscErrorCode dmplex_get_face_connectivity(DM plex, PetscInt nfaces, PetscInt sloc, PetscInt *sndglno);
