#include "DMPlex_Reader_C.h"

PetscErrorCode dmplex_mark_halo_regions(DM *plex)
{
  DMLabel         lblDepth, lblHalo;
  PetscInt        p, pStart, pEnd, cStart, cEnd, nleaves, nroots, npoints;
  PetscInt        a, adjSize, *adj = NULL, cl, clSize, *closure = NULL;
  const PetscInt *ilocal, *points;
  PetscSF         pointSF;
  IS              haloPoints;
  PetscBool       hasValue, useCone;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetChart(*plex, &pStart, &pEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(*plex, 0,  &cStart, &cEnd);CHKERRQ(ierr);
  /* Use star(closure(p)) adjacency to get neighbouring cells */
  ierr = DMPlexGetAdjacencyUseCone(*plex, &useCone);CHKERRQ(ierr);
  ierr = DMPlexSetAdjacencyUseCone(*plex, PETSC_TRUE);CHKERRQ(ierr);

  /* Create labels for L1 and L2 halos */
  ierr = DMPlexGetLabel(*plex, "depth", &lblDepth);CHKERRQ(ierr);
  ierr = DMPlexCreateLabel(*plex, "HaloRegions");CHKERRQ(ierr);
  ierr = DMPlexGetLabel(*plex, "HaloRegions", &lblHalo);CHKERRQ(ierr);
  ierr = DMLabelCreateIndex(lblHalo, pStart, pEnd);CHKERRQ(ierr);

  /* Loop over point SF and mark the entire halo region L2 */
  ierr = DMGetPointSF(*plex, &pointSF);CHKERRQ(ierr);
  ierr = PetscSFGetGraph(pointSF, &nroots, &nleaves, &ilocal, NULL);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {ierr = DMLabelSetValue(lblHalo, ilocal[p], 2);CHKERRQ(ierr);}

  /* Loop over halo cells and test for non-halo neighbouring cells */
  ierr = DMLabelGetStratumIS(lblHalo, 2, &haloPoints);CHKERRQ(ierr);
  ierr = DMLabelGetStratumSize(lblHalo, 2, &npoints);CHKERRQ(ierr);
  ierr = ISGetIndices(haloPoints, &points);CHKERRQ(ierr);
  for (p = 0; p < npoints; p++) {
    const PetscInt cell = points[p];
    if (cStart <= cell && cell < cEnd) {
      adjSize = PETSC_DETERMINE;
      ierr = DMPlexGetAdjacency(*plex, cell, &adjSize, &adj);CHKERRQ(ierr);
      for (a = 0; a < adjSize; ++a) {
        const PetscInt neigh = adj[a];
        if (neigh != cell && cStart <= neigh && neigh < cEnd) {
          /* If the neighbouring cell is not in L2; mark this cell as L1 */
          ierr = DMLabelStratumHasPoint(lblHalo, 2, neigh, &hasValue);CHKERRQ(ierr);
          if (!hasValue) {
            ierr = DMPlexGetTransitiveClosure(*plex, cell, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
            for (cl = 0; cl < 2*clSize; cl+=2) {
              /* L1 is a subset of L2 */
              ierr = DMLabelStratumHasPoint(lblHalo, 2, closure[cl], &hasValue);CHKERRQ(ierr);
              if (hasValue) {ierr = DMLabelSetValue(lblHalo, closure[cl], 1);CHKERRQ(ierr);}
            }
          }
        }
      }
    }
  }
  ierr = ISRestoreIndices(haloPoints, &points);CHKERRQ(ierr);
  ierr = ISDestroy(&haloPoints);CHKERRQ(ierr);
  if (adj) {ierr = PetscFree(adj);CHKERRQ(ierr);}
  if (closure) {ierr = DMPlexRestoreTransitiveClosure(*plex, 0, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);}

  ierr = DMPlexSetAdjacencyUseCone(*plex, useCone);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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
