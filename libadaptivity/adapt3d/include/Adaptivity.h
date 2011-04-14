/*
  Copyright (C) 2006 Imperial College London and others.

  Please see the AUTHORS file in the main source directory for a full list
  of copyright holders.

  Gerard Gorman
  Applied Modelling and Computation Group
  Department of Earth Science and Engineering
  Imperial College London

  adrian@Imperial.ac.uk

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  USA
*/

#ifndef ADAPTIVITY_H
#define ADAPTIVITY_H

#include "confdefs.h"

#ifdef USING_DOUBLE_PRECISION
typedef double afloat_t;
#else
typedef float afloat_t;
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <deque>
#include <map>
#include <string>
#include <vector>
#include <set>

#include <string.h>

#include "vtk.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/// C++ class wrapper for serial adaptivity
class Adaptivity{
 public:
  /// Constructor
  Adaptivity();

  /// What percentage of edges have length outside the interval
  /// [1-w, 1+w]? This can be used to decide whether to adapt
  /// or not.
  afloat_t edgeLengthDistribution(afloat_t w);

  /// Perform adapt
  void adapt();

  /// Enable/disable automatic detection and protection of co-planar patches.
  void enableGeometryDiscovery();

  /// Enable/disable node movement (smoothening).
  void enableMovementConnectivity();

  /// Enable/disable mesh refinement
  void enableRefinement();

  /// Enable/disable surface locking
  void enableSurfaceLock();

  /// Enable/disable automatic detection and protection of co-planar patches.
  void disableGeometryDiscovery();

  /// Enable/disable node movement (smoothening).
  void disableMovementConnectivity();

  /// Enable/disable mesh refinement;
  void disableRefinement();

  /// Enable/disable surface locking
  void disableSurfaceLock();

#ifdef HAVE_VTK
  /// Get a vtkUnstructuredGrid object that contains the adapted mesh.
  vtkUnstructuredGrid* get_adapted_vtu();
#endif

  /// Get mesh partition halo (parallel only).
  void getHalo(int *_numPrivateNodesint, int *_ATOSEN, int *_ATOREC, int *_Gather, int *_Scatter);

  /// Get the gimension of the mesh
  void getMeshDimensions(int *_NNodes, int *_NElements, int *_NSElements);
  
  /// Get the surface lables.
  void get_surface_ids(std::vector<int> &sids);

  /// Get the surface mesh.
  void get_surface_mesh(std::vector<int> &SENList);

  /// Get the points of the mesh.
  void getVertices(afloat_t *_X, afloat_t *_Y, afloat_t *Z);

  /// Get the element lables.
  void getVolumeID(int *_newvolumeID);

  /// Get the volume mesh
  void getVolumeMesh(int *_newENList);

  /// Set functional tolerances.
  void setFunctionalTolerance(afloat_t);

  /// Set the partition halo (parallel only).
  void setHalo(int _numPrivateNodesint, int *_ATOSEN, int *_ATOREC, int *_Gather, int *_Scatter);

  /// Set the maximum number of adaptive sweeps through the mesh.
  void set_adapt_sweeps(int _iterations);

  /// Set the metic tensor field.
  void set_metric(afloat_t *_Metric);

  /// Set the surface labels/ids.
  void set_surface_ids(std::vector<int> &sids);

  /// Set the surface mesh.
  void set_surface_mesh(std::vector<int> &SENList);

  /// Set the mesh points.
  void set_points(afloat_t *X, afloat_t *Y, afloat_t *Z, int _NNodes);

  /// Set the element lables.
  void setVolumeID(int *_volumeID);

  /// Set the volume mesh.
  void set_volume_mesh(int *_ENList,   int _NElements);

#ifdef HAVE_VTK
  /// Set the input from a vtkUnstructuredGrid. Set interpolate=true
  /// to interpolate all the fields.
  void set_from_vtk(vtkUnstructuredGrid *ug, bool interpolate);
#endif

  /// Turn on/off verbose mode.
  void verbose_off();

  /// Turn on/off verbose mode.
  void verbose_on();
 private:
  static bool verbose;

  int getNProcessors() const;
  double volume(int n1, int n2, int n3, int n4) const;
  double volume(int n1, int n2, int n3, afloat_t x, afloat_t y, afloat_t z) const;
  double volume(afloat_t x1, afloat_t y1, afloat_t z1,
		afloat_t x2, afloat_t y2, afloat_t z2,
		afloat_t x3, afloat_t y3, afloat_t z3,
		afloat_t x4, afloat_t y4, afloat_t z4) const;

  int newNNodes, newNElements, newNSElements;

#ifdef HAVE_VTK  
  vtkUnstructuredGrid *ug;
#else
  void *ug;
#endif

  // Initial mesh.
  std::vector<afloat_t> X, Y, Z;
  std::vector<int> ENList, SENList, surfID;
  const int *volumeID;
  int numberingOffset;
  
  //New mesh.
  afloat_t *newX, *newY, *newZ;
  
  // Extra stuff required for parallel computations
  int NPrivateNodes, NProcs;
  const int *ATOSEN, *ATOREC;
  const int *Gather, *Scatter;
  int NGather, NScatter;

  int newNPrivateNodes;
  std::vector<int> newATOSEN, newATOREC;
  std::vector<int> newGather, newScatter;

  // Metric being used for directing mesh adaptivity.
  std::vector<afloat_t> Metric;
  
  // Adaptivity options
  int MaxNumberAdaptIterations;
  int AdaptOpts[6];

  // even more adaptivity options
  // arg        | default | comment
  //---------------------------------------------------------------------------------------
  int SRFGMY; // .FALSE. | is TRUE if the surface mesh should be kept intact during adapting
  int CLCGMY; // .TRUE.  | is TRUE if the geometry should be calculated, and ignore SNLIST
  int TWOSTG; // .FALSE. | two stages of adapting, with no refinement on first
  int TOGTHR; // .TRUE.  | lumps node movement adaptivity in with connectivity changes
  int chcnsy; // .FALSE. | do consistancy checks
  int dbg;    // .FALSE. | debug adaptivity
  
  int nloc, snloc;

  // MESTP1 is the value of functional after which elements are
  // considered for adaption changes. 0.5 to 1.0 seems to work ok. The
  // value can be negative.  0.5 corresponds to a minimum insphere
  // radius 0f 0.3 relative 1 and max-edge size 2.
  afloat_t MESTP1;

  // The new mesh will be stored in these buffers.
  std::vector<int> intBuffer;
  std::vector<afloat_t> floatBuffer;
  
  // Store fields.
  int nfields;
  std::vector<afloat_t> fields;
  std::vector<int> nfreedom;
    
  // Interpolate fields
  bool interpolate;

  // These base pointers store the new mesh in the buffers
  int 
    NWSZEN, NWSZSN, NWSZNN, NWNDLC, NWSROW,
    NWENLB, NWENLS, NWSNLB, NWSNLS, NWSFID,
    NWELRG, NWNODX, NWNODY, NWNODZ,
    NEWMTX, NEWFLD,
    ADPBIG, ADPNOD;

  int UsingSurfaces;
};

#endif
