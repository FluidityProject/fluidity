#include "DMPlex_Reader_C.h"

PetscErrorCode dmplex_get_mesh_connectivity(DM plex, PetscInt nnodes, PetscInt loc, PetscInt *ndglno)
{
  PetscInt c, cStart, cEnd, vStart, vEnd, idx, ci, nclosure, point;
  PetscInt *closure=NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetHeightStratum(plex, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  for (idx=0, c=cStart; c<cEnd; c++) {
    ierr = DMPlexGetTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
    for (ci=0; ci<nclosure; ci++) {
      point = closure[2*ci];
      if (vStart <= point && point < vEnd) ndglno[idx++] = point - vStart + 1;
    }
  }
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_face_connectivity(DM plex, PetscInt nfaces, PetscInt sloc, PetscInt *sndglno)
{
  PetscInt f, fStart, fEnd, vStart, vEnd, idx, ci, nclosure, point;
  PetscInt *closure=NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetHeightStratum(plex, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  for (idx=0, f=fStart; f<fEnd; f++) {
    ierr = DMPlexGetTransitiveClosure(plex, f, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
    for (ci=0; ci<nclosure; ci++) {
      point = closure[2*ci];
      if (vStart <= point && point < vEnd) sndglno[idx++] = point - vStart + 1;
    }
  }
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, f, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
