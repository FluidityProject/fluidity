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

#ifndef DISCRETEGEOMETRYCONSTRAINTS_H
#define DISCRETEGEOMETRYCONSTRAINTS_H

#include "confdefs.h"

#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <string.h>

#include "vtk.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/** This class is used to: identify the boundary of the domain;
    uniquely label connected co-linear patches of surface elements
    (these can be used to prevent adaptivity coarsening these patches
    and smoothening out features); evaluate a characteristic length
    scale for these patches (these "constraints" can be added to the
    metric tensor field before gradation is applied in order to get
    good quality elements near the geometry).
*/

class DiscreteGeometryConstraints{
 public:
  /// Default constructor.
  DiscreteGeometryConstraints();

  /// Get constraints.
  /// @param[out] Maximum lengths - 3x3 symmetric tensor field.
  void get_constraints(std::vector<double> &max_len);

  /// Get copy of coplanar_ids.
  /// @param[out] Surface element co-planar patch id's.
  void get_coplanar_ids(std::vector<int> &ids);

  /// Get surface element-node list.
  /// @param[out] ENList element-node list
  void get_surface(std::vector<int> &ENList);

  /// Set dot product tolerence - used to decide if elements are co-planar
  void set_coplanar_tolerance(double tol);

  /// Set vtkUnstructuredGrid as the input mesh.
  void set_volume_input(vtkUnstructuredGrid *ug);

  /// Set vtkUnstructuredGrid as the input mesh.
  void set_surface_input(vtkUnstructuredGrid *ug, std::vector<int> &ENList, std::vector<int> &ids);

  /// Turn off verbose mode.
  void verbose_off();

  /// Turn on verbose mode.
  void verbose_on();

#ifdef HAVE_VTK
  /// Write out VTU of the constraints.
  void write_vtk(std::string filename);
#endif

 private:
  /// Calculate co-planar patches.
  void get_coplanar_ids();

  /// Discover limitations on maximum edge length at the surface of
  /// the geometry. This is a simplification of algorithm in Fluidity
  /// and as such has a few short comings. The full method will be implemented later.
  void find_geometry_constraints();

  /// Normal to a triangle.
  void normal(int v1, int v2, int v3, double *n);

  static bool verbose;
  vtkUnstructuredGrid *ug;
  std::vector<int> SENList, EEList, coplanar_ids;
  double COPLANAR_MAGIC_NUMBER;
  std::map< int, std::set<size_t> > SNEList;
  std::vector<double> gconstraint;
};

#endif
