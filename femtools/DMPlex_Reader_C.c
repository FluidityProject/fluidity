#include "DMPlex_Reader_C.h"

PetscErrorCode dmplex_derive_loc(DM *plex, PetscInt *point, PetscInt *loc)
{
  PetscInt vStart, vEnd, cl, clSize, *closure = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(*plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetTransitiveClosure(*plex, *point, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
  for (*loc = 0, cl = 0; cl < 2*clSize; cl+=2) {if (vStart <= closure[cl] && closure[cl] < vEnd) (*loc)++;}
  ierr = DMPlexRestoreTransitiveClosure(*plex, *point, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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

PetscErrorCode dmplex_get_reordering(DM *plex, PetscInt depth, IS *permutation, IS *reordering)
{
  MPI_Comm        comm;
  PetscInt        v, p, pStart, pEnd, npoints, idx, size, *ordering;
  DMLabel         lblHalo;
  IS              haloL1, haloL2;
  const PetscInt *points, *perm;
  PetscBool       hasPoint;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) *plex, &comm);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(*plex, "HaloRegions", &lblHalo);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(*plex, depth, &pStart, &pEnd);CHKERRQ(ierr);

  ierr = PetscMalloc1(pEnd - pStart, &ordering);CHKERRQ(ierr);
  ierr = ISGetLocalSize(*permutation, &size);CHKERRQ(ierr);
  ierr = ISGetIndices(*permutation, &perm);CHKERRQ(ierr);

  if (!lblHalo) {
    for (idx = 0, v = 0; v < size; v++) {
      const PetscInt point = perm[v];
      if (pStart <= point && point < pEnd) ordering[point - pStart] = idx++;
    }
    ierr = ISRestoreIndices(*permutation, &perm);CHKERRQ(ierr);
    ierr = ISCreateGeneral(comm, pEnd-pStart, ordering, PETSC_OWN_POINTER, reordering);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /* Add owned points first */
  ierr = PetscMalloc1(pEnd - pStart, &ordering);CHKERRQ(ierr);
  for (idx = 0, v = 0; v < size; v++) {
    const PetscInt point = perm[v];
    if (pStart <= point && point < pEnd) {
      ierr = DMLabelHasPoint(lblHalo, point, &hasPoint);CHKERRQ(ierr);
      if (!hasPoint) ordering[point - pStart] = idx++;
    }
  }

  /* Add entities in L1 halo region */
  ierr = DMLabelGetStratumIS(lblHalo, 1, &haloL1);CHKERRQ(ierr);
  ierr = DMLabelGetStratumSize(lblHalo, 1, &npoints);CHKERRQ(ierr);
  ierr = ISGetIndices(haloL1, &points);CHKERRQ(ierr);
  for (v = 0; v < size; v++) {
    const PetscInt point = perm[v];
    if (pStart <= point && point < pEnd) {
      ierr = DMLabelStratumHasPoint(lblHalo, 1, point, &hasPoint);CHKERRQ(ierr);
      if (hasPoint) ordering[point - pStart] = idx++;
    }
  }
  ierr = ISRestoreIndices(haloL1, &points);CHKERRQ(ierr);
  ierr = ISDestroy(&haloL1);CHKERRQ(ierr);

  /* Add entities in L2 halo region, but not already in L1 */
  ierr = DMLabelGetStratumIS(lblHalo, 2, &haloL2);CHKERRQ(ierr);
  ierr = DMLabelGetStratumSize(lblHalo, 2, &npoints);CHKERRQ(ierr);
  for (v = 0; v < size; v++) {
    const PetscInt point = perm[v];
    if (pStart <= point && point < pEnd) {
      ierr = DMLabelStratumHasPoint(lblHalo, 1, point, &hasPoint);CHKERRQ(ierr);
      if (hasPoint) continue;
      ierr = DMLabelStratumHasPoint(lblHalo, 2, point, &hasPoint);CHKERRQ(ierr);
      if (hasPoint) ordering[point - pStart] = idx++;
    }
  }
  ierr = ISDestroy(&haloL2);CHKERRQ(ierr);
  ierr = ISRestoreIndices(*permutation, &perm);CHKERRQ(ierr);

  ierr = ISCreateGeneral(comm, pEnd-pStart, ordering, PETSC_OWN_POINTER, reordering);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_mesh_connectivity(DM plex, PetscInt nnodes, PetscInt loc,
                                            IS *rnbrCells, IS *rnbrVertices, PetscInt *ndglno)
{
  PetscInt c, cStart, cEnd, vStart, vEnd, idx, ci, nclosure, point;
  PetscInt *closure=NULL;
  const PetscInt *cells, *vertices;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetHeightStratum(plex, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = ISGetIndices(*rnbrCells, &cells);CHKERRQ(ierr);
  ierr = ISGetIndices(*rnbrVertices, &vertices);CHKERRQ(ierr);
  for (c=cStart; c<cEnd; c++) {
    ierr = DMPlexGetTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
    for (idx=0, ci=0; ci<nclosure; ci++) {
      point = closure[2*ci];
      if (vStart <= point && point < vEnd) ndglno[loc*cells[c] + idx++] = vertices[point - vStart] + 1;
    }
  }
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(*rnbrCells, &cells);CHKERRQ(ierr);
  ierr = ISRestoreIndices(*rnbrVertices, &vertices);CHKERRQ(ierr);
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
  if (!label) {*nfacets = 0; PetscFunctionReturn(0);}
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

PetscErrorCode dmplex_get_surface_connectivity(DM plex, const char label_name[], PetscInt nfacets, PetscInt sloc, IS *rnbrVertices, PetscInt *sndglno, PetscInt *boundary_ids)
{
  PetscInt        v, vStart, vEnd, nvalues, p, npoints, idx, ci, nclosure, vertex, nvertices;
  const PetscInt *values, *points, *vertices;
  DMLabel         label;
  IS              valueIS, pointIS;
  PetscInt       *closure = NULL;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(plex, label_name, &label);CHKERRQ(ierr);
  if (!label) PetscFunctionReturn(0);
  ierr = DMLabelGetNumValues(label, &nvalues);CHKERRQ(ierr);
  ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
  ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
  ierr = ISGetIndices(*rnbrVertices, &vertices);CHKERRQ(ierr);
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
          sndglno[idx*sloc+nvertices++] = vertices[vertex - vStart] + 1;
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
    ierr = DMPlexRestoreTransitiveClosure(plex, 0, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(*rnbrVertices, &vertices);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_halo_receives(DM *plex, PetscInt nprocs, PetscInt height, IS *renumbering,
                                        PetscInt *nrecv_l1, IS *receives_l1,
                                        PetscInt *nrecv_l2, IS *receives_l2) {
  MPI_Comm           comm;
  PetscSF            sfPoint;
  DMLabel            lblHalo;
  PetscInt           p, pStart, pEnd, nleaves, size_l1, size_l2, offset;
  PetscInt          *recv_arr_l1, *recv_arr_l2, *idx_l1, *idx_l2;
  PetscSection       recvSectionL1, recvSectionL2;
  const PetscInt    *ilocal, *perm;
  const PetscSFNode *iremote;
  PetscBool          hasPointL1, hasPointL2;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) *plex, &comm);CHKERRQ(ierr);
  ierr = DMGetPointSF(*plex, &sfPoint);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(*plex, "HaloRegions", &lblHalo);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(*plex, height, &pStart, &pEnd);CHKERRQ(ierr);
  ierr = ISGetIndices(*renumbering, &perm);CHKERRQ(ierr);

  /* Count L1/L2 receives using sections */
  ierr = PetscSFGetGraph(sfPoint, NULL, &nleaves, &ilocal, &iremote);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &recvSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(recvSectionL1, 0, nprocs);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &recvSectionL2);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(recvSectionL2, 0, nprocs);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {
    if (pStart <= ilocal[p] && ilocal[p] < pEnd) {
      const PetscInt proc = iremote[p].rank;
      ierr = DMLabelStratumHasPoint(lblHalo, 1, ilocal[p], &hasPointL1);CHKERRQ(ierr);
      if(hasPointL1) {ierr = PetscSectionAddDof(recvSectionL1, proc, 1);CHKERRQ(ierr);}
      ierr = DMLabelStratumHasPoint(lblHalo, 2, ilocal[p], &hasPointL2);CHKERRQ(ierr);
      if(hasPointL2) {ierr = PetscSectionAddDof(recvSectionL2, proc, 1);CHKERRQ(ierr);}
    }
  }
  /* Set the number of receives per process */
  ierr = PetscSectionSetUp(recvSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionSetUp(recvSectionL2);CHKERRQ(ierr);
  for (p = 0; p < nprocs; p++) {
    ierr = PetscSectionGetDof(recvSectionL1, p, &(nrecv_l1[p]));CHKERRQ(ierr);
    ierr = PetscSectionGetDof(recvSectionL2, p, &(nrecv_l2[p]));CHKERRQ(ierr);
  }
  /* Allocate and fill the receive arrays */
  ierr = PetscSectionGetStorageSize(recvSectionL1, &size_l1);CHKERRQ(ierr);
  ierr = PetscMalloc1(size_l1, &recv_arr_l1);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(recvSectionL2, &size_l2);CHKERRQ(ierr);
  ierr = PetscMalloc1(size_l2, &recv_arr_l2);CHKERRQ(ierr);
  ierr = PetscCalloc2(nprocs, &idx_l1, nprocs, &idx_l2);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {
    if (pStart <= ilocal[p] && ilocal[p] < pEnd) {
      const PetscInt proc = iremote[p].rank;
      const PetscInt point = perm[ilocal[p] - pStart]; 
      ierr = DMLabelStratumHasPoint(lblHalo, 1, ilocal[p], &hasPointL1);CHKERRQ(ierr);
      ierr = DMLabelStratumHasPoint(lblHalo, 2, ilocal[p], &hasPointL2);CHKERRQ(ierr);
      if(hasPointL1) {
        /* Add L1 points to both receive lists */
        ierr = PetscSectionGetOffset(recvSectionL1, proc, &offset);CHKERRQ(ierr);
        recv_arr_l1[offset+idx_l1[proc]] = point;
        ierr = PetscSectionGetOffset(recvSectionL2, proc, &offset);CHKERRQ(ierr);
        recv_arr_l2[offset+idx_l1[proc]] = point;
        idx_l1[proc]++;
      } else if (hasPointL2) {
        /* Add L2 points at the end of the proc-specific segment of the receive array */
        ierr = PetscSectionGetOffset(recvSectionL2, proc, &offset);CHKERRQ(ierr);
        recv_arr_l2[offset+nrecv_l1[proc]+idx_l2[proc]] = point;
        idx_l2[proc]++;
      }
    }
  }
  /* Return receive arrays wrapped as IS objects */
  ierr = ISCreateGeneral(comm, size_l1, recv_arr_l1, PETSC_OWN_POINTER, receives_l1);CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm, size_l2, recv_arr_l2, PETSC_OWN_POINTER, receives_l2);CHKERRQ(ierr);

  ierr = ISRestoreIndices(*renumbering, &perm);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&recvSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&recvSectionL2);CHKERRQ(ierr);
  ierr = PetscFree2(idx_l1, idx_l2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_halo_sends(DM *plex, PetscInt nprocs, PetscInt height, IS *renumbering,
                                     PetscInt *nsend_l1, IS *sends_l1,
                                     PetscInt *nsend_l2, IS *sends_l2) {
  MPI_Comm           comm;
  PetscSF            sfPoint;
  DMLabel            lblHalo;
  PetscInt           p, pStart, pEnd, size, size_l1, size_l2, roff, offset, r, nranks;
  PetscInt          *leafHaloLevel, *rootHaloLevel, *send_arr_l1, *send_arr_l2, *idx_l1, *idx_l2;
  PetscSection       rootSection, leafSection, sendSectionL1, sendSectionL2;
  IS                 rootRanks, leafRanks;
  const PetscInt    *perm, *rranks;
  PetscBool          hasPointL1, hasPointL2;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) *plex, &comm);CHKERRQ(ierr);
  ierr = DMGetPointSF(*plex, &sfPoint);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(*plex, "HaloRegions", &lblHalo);CHKERRQ(ierr);
  /* ierr = DMPlexGetDepthStratum(*plex, height, &pStart, &pEnd);CHKERRQ(ierr); */
  ierr = ISGetIndices(*renumbering, &perm);CHKERRQ(ierr);

  /* Determine send targets for each root */
  ierr = PetscSectionCreate(comm, &rootSection);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &leafSection);CHKERRQ(ierr);
  ierr = DMPlexDistributeOwnership(*plex, rootSection, &rootRanks, leafSection, &leafRanks);CHKERRQ(ierr);
  ierr = ISGetIndices(rootRanks, &rranks);CHKERRQ(ierr);
  /* In order to create the appropriate send list ordering we need to know the halo level
     of the target node. To derive this, we PetscSFGather the halo level over the point SF. */
  ierr = DMPlexGetChart(*plex, &pStart, &pEnd);CHKERRQ(ierr);
  ierr = PetscCalloc1(pEnd-pStart, &leafHaloLevel);CHKERRQ(ierr);
  PetscInt nleaves;
  const PetscInt    *ilocal;
  const PetscSFNode *iremote;
  ierr = PetscSFGetGraph(sfPoint, NULL, &nleaves, &ilocal, &iremote);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {
    /* PetscSFGather ignore the local point indirection, so we have to apply it manually */
    ierr = DMLabelStratumHasPoint(lblHalo, 1, ilocal[p], &hasPointL1);CHKERRQ(ierr);
    ierr = DMLabelStratumHasPoint(lblHalo, 2, ilocal[p], &hasPointL2);CHKERRQ(ierr);
    if(hasPointL1) leafHaloLevel[p] = 1;
    else if(hasPointL2) leafHaloLevel[p] = 2;
  }
  ierr = PetscSectionGetStorageSize(rootSection, &size);CHKERRQ(ierr);
  ierr = PetscMalloc1(size, &rootHaloLevel);CHKERRQ(ierr);
  ierr = PetscSFGatherBegin(sfPoint, MPIU_INT, leafHaloLevel, rootHaloLevel);CHKERRQ(ierr);
  ierr = PetscSFGatherEnd(sfPoint, MPIU_INT, leafHaloLevel, rootHaloLevel);CHKERRQ(ierr);

  /* Now we're only interested in one stratum */
  ierr = DMPlexGetDepthStratum(*plex, height, &pStart, &pEnd);CHKERRQ(ierr);

  /* Count L1/L2 sends using sections */
  ierr = PetscSectionCreate(comm, &sendSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(sendSectionL1, 0, nprocs);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &sendSectionL2);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(sendSectionL2, 0, nprocs);CHKERRQ(ierr);
  for (p = pStart; p < pEnd; p++) {
    ierr = PetscSectionGetOffset(rootSection, p, &offset);CHKERRQ(ierr);
    ierr = PetscSectionGetDof(rootSection, p, &nranks);CHKERRQ(ierr);
    for (r = 0; r < nranks; r++) {
      const PetscInt proc = rranks[offset+r];
      if (rootHaloLevel[offset+r] == 1) {
        ierr = PetscSectionAddDof(sendSectionL1, proc, 1);CHKERRQ(ierr);
        ierr = PetscSectionAddDof(sendSectionL2, proc, 1);CHKERRQ(ierr);
      } else if (rootHaloLevel[offset+r] == 2) {
        ierr = PetscSectionAddDof(sendSectionL2, proc, 1);CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscSectionSetUp(sendSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionSetUp(sendSectionL2);CHKERRQ(ierr);

  for (p = 0; p < nprocs; p++) {
    ierr = PetscSectionGetDof(sendSectionL1, p, &(nsend_l1[p]));CHKERRQ(ierr);
    ierr = PetscSectionGetDof(sendSectionL2, p, &(nsend_l2[p]));CHKERRQ(ierr);
  }
  /* Allocate and fill the send arrays */
  ierr = PetscSectionGetStorageSize(sendSectionL1, &size_l1);CHKERRQ(ierr);
  ierr = PetscMalloc1(size_l1, &send_arr_l1);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(sendSectionL2, &size_l2);CHKERRQ(ierr);
  ierr = PetscMalloc1(size_l2, &send_arr_l2);CHKERRQ(ierr);
  ierr = PetscCalloc2(nprocs, &idx_l1, nprocs, &idx_l2);CHKERRQ(ierr);
  for (p = pStart; p < pEnd; p++) {
    ierr = PetscSectionGetOffset(rootSection, p, &roff);CHKERRQ(ierr);
    ierr = PetscSectionGetDof(rootSection, p, &nranks);CHKERRQ(ierr);
    for (r = 0; r < nranks; r++) {
      const PetscInt proc = rranks[roff+r];
      const PetscInt haloLevel = rootHaloLevel[roff+r];
      const PetscInt point = perm[p - pStart];
      if (haloLevel == 1) {
        ierr = PetscSectionGetOffset(sendSectionL1, proc, &offset);CHKERRQ(ierr);
        send_arr_l1[offset+idx_l1[proc]] = point;
        ierr = PetscSectionGetOffset(sendSectionL2, proc, &offset);CHKERRQ(ierr);
        send_arr_l2[offset+idx_l1[proc]] = point;
        idx_l1[proc]++;
      } else if (haloLevel == 2) {
        ierr = PetscSectionGetOffset(sendSectionL2, proc, &offset);CHKERRQ(ierr);
        send_arr_l2[offset+nsend_l1[proc]+idx_l2[proc]] = point;
        idx_l2[proc]++;
      }
    }
  }
  /* Return send arrays wrapped as IS objects */
  ierr = ISCreateGeneral(comm, size_l1, send_arr_l1, PETSC_OWN_POINTER, sends_l1);CHKERRQ(ierr);
  ierr = ISCreateGeneral(comm, size_l2, send_arr_l2, PETSC_OWN_POINTER, sends_l2);CHKERRQ(ierr);

  ierr = PetscSectionDestroy(&rootSection);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&leafSection);CHKERRQ(ierr);
  ierr = ISDestroy(&rootRanks);CHKERRQ(ierr);
  ierr = ISDestroy(&leafRanks);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&sendSectionL1);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&sendSectionL1);CHKERRQ(ierr);
  ierr = PetscFree(leafHaloLevel);CHKERRQ(ierr);
  ierr = PetscFree2(idx_l1, idx_l2);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
