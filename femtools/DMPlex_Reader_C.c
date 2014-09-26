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

PetscErrorCode dmplex_get_num_surface_facets(DM plex, const char label_name[], PetscInt *nfacets)
{
  PetscInt v, nvalues, npoints;
  const PetscInt *values;
  DMLabel label;
  IS valueIS;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetLabel(plex, label_name, &label);CHKERRQ(ierr);
  ierr = DMLabelGetNumValues(label, &nvalues);CHKERRQ(ierr);
  ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
  ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
  for (*nfacets=0, v=0; v<nvalues; v++) {
    ierr = DMLabelGetStratumSize(label, values[v], &npoints);CHKERRQ(ierr);
    *nfacets += npoints;
  }
  ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
  ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_surface_connectivity(DM plex, const char label_name[], PetscInt nfacets, PetscInt sloc, PetscInt *sndglno, PetscInt *boundary_ids)
{
  PetscInt        v, vStart, vEnd, nvalues, p, npoints, idx, ci, nclosure, vertex, nvertices;
  const PetscInt *values, *points;
  DMLabel         label;
  IS              valueIS, pointIS;
  PetscInt       *closure = NULL;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(plex, label_name, &label);CHKERRQ(ierr);
  ierr = DMLabelGetNumValues(label, &nvalues);CHKERRQ(ierr);
  ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
  ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
  /* Loop over all marker values in the supplied label */
  for (idx = 0 , v = 0; v < nvalues; v++) {
    ierr = DMLabelGetStratumSize(label, values[v], &npoints);CHKERRQ(ierr);
    ierr = DMLabelGetStratumIS(label, values[v], &pointIS);CHKERRQ(ierr);
    ierr = ISGetIndices(pointIS, &points);CHKERRQ(ierr);
    for (p = 0; p < npoints; p++) {
      /* Derive vertices for each marked facet */
      ierr = DMPlexGetTransitiveClosure(plex, points[p], PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
      for (nvertices = 0, ci = 0; ci < nclosure; ci++) {
        vertex = closure[2*ci];
        if (vStart <= vertex && vertex < vEnd) {
          sndglno[idx*sloc+nvertices++] = vertex - vStart + 1;
        }
      }
      /* Store associated boundary ID */
      boundary_ids[idx] = values[v];
      idx++;
    }
    ierr = ISRestoreIndices(pointIS, &points);CHKERRQ(ierr);
    ierr = ISDestroy(&pointIS);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
  ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, points[p], PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
