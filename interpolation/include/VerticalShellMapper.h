/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
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
#ifndef VERTICALSHELLMAPPER_H
#define VERTICALSHELLMAPPER_H

#include "confdefs.h"

#include "projections.h"
#include "spatialindex/SpatialIndex.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "vtk.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>


class VerticalShellMapper{
 public:
  VerticalShellMapper();
  ~VerticalShellMapper();
  
  void DumpVTU(const char *) const;
  int Find(double x, double y, double z, double *shape);
  int pFind(const double *x, const double *y, const double *z, int cnt,
            int *tid, double *shape, int *rank);
  bool IsParallel() const;
  void SetShell(const double *, const double *, const double *, const int *, int, int);
  void SetSphericalGeometry(bool);
  
 private:
  bool verbose;
  size_t data_block_size;
  bool spherical_geometry;
  std::vector<double> sx, sy;
  std::vector<int>
    Triangles_north, Triangle_ids_north,
    Triangles_south, Triangle_ids_south,
    Triangles;
  
  // I'm going to guess that only one storage manager is necessary -
  // and data is seperated using the indexId. Get it working for 2
  // first and see then if it can be reduced.
  SpatialIndex::IStorageManager *storage_manager_north, *storage_manager_south, *storage_manager;
  SpatialIndex::ISpatialIndex *rtree_north, *rtree_south, *rtree;
  SpatialIndex::id_type indexId_north, indexId_south, indexId;
};

#endif
