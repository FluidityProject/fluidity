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
#ifndef PICKER_H
#define PICKER_H

#include "confdefs.h"

#include "projections.h"
#include "spatialindex/SpatialIndex.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vtk.h>

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>

// Customised version of MyDataStream class in
// regressiontest/rtree/RTreeBulkLoad.cc in spatialindex 1.4.0
class PickerMeshDataStream : public SpatialIndex::IDataStream{
 public:
  PickerMeshDataStream(const double *x, const double *y, const double *z,
                       const int *enlist, int zero_index, int nelements, int loc);

  PickerMeshDataStream(const double *x, const double *y, const double *z,
                       const int *enlist, int zero_index, int nelements, int loc,
                       const int *live_set);

  PickerMeshDataStream(const double *x, const double *y,
                       const int *enlist, int zero_index, int nelements, int loc);

  PickerMeshDataStream(const double *x, const double *y,
                       const int *enlist, int zero_index, int nelements, int loc,
                       const int *live_set);

  virtual ~PickerMeshDataStream();
  
  virtual SpatialIndex::IData* getNext();
  void readNextEntry();

  virtual bool hasNext();

#ifdef EXPERIMENTAL_SPATIALINDEX
  virtual uint32_t size();
#else
  virtual size_t size();
#endif
  virtual void rewind();
 protected:
  SpatialIndex::RTree::Data* m_pNext;
  void init(const double *x, const double *y, const double *z,
            const int *enlist, int zero_index, int nelements, int loc,
            const int *live_set);
  const double *x, *y, *z, *pos[3];
  const int *enlist, *live_set;
  int dim, index, nelements, loc, zero_index;
  
  size_t data_block_size;
};

class Picker{
 public:
  Picker(const double& tol = 1.0e-3);
  ~Picker();

  void DumpVTU(const char *) const;
  
  void SetOwnershipTolerance(const double& tol);
  double GetOwnershipTolerance() const;
  int Find(const double *x, int *eid, double *shape);
  int Find(double x, double y, int *eid, double *shape);
  int Find(double x, double y, double z, int *eid, double *shape);

  int pFind(const double *x, const double *y, const double *z, int cnt,
            int *tid, double *shape, int *rank, int different_domains);
  
  bool IsParallel() const;

  // Set z = NULL if no z coordinate -- ie if 2D
  void SetSource(const double *x, const double *y,
                 const int *elements, int nelements,
                 int _element_type, int zero_index);

  void SetSource(const double *x, const double *y, const double *z,
                 const int *enlist, int nelements,
                 int _element_type, int zero_index);
  
  void SetSphericalGeometry(bool);
  
  void VerboseOff();
  void VerboseOn();
 private:
  void CreateStorageManagers();

  bool verbose;
  unsigned long dimension, edimension, element_type, nloc;
  double tol;
  size_t data_block_size;
  bool spherical_geometry;
  std::vector<double> sx, sy, sz;
  const double *_x, *_y, *_z;
  
  std::vector<int>
    elements_north, element_ids_north,
    elements_south, element_ids_south,
    elements;
  
  // I'm going to guess that only one storage manager is necessary -
  // and data is seperated using the indexId. Get it working for 2
  // first and see then if it can be reduced.
  SpatialIndex::IStorageManager *storage_manager_north, *storage_manager_south, *storage_manager;
  SpatialIndex::ISpatialIndex *rtree_north, *rtree_south, *rtree;
  SpatialIndex::id_type indexId_north, indexId_south, indexId;

  const unsigned long TRI, QUAD, TETRA, HEX;

  // R-Tree parameters - need to see how these need to be optimized
  // for 2D and 3D.
  double fillFactor;
  unsigned long indexCapacity;
  unsigned long leafCapacity;
  SpatialIndex::RTree::RTreeVariant variant;
};

#endif
