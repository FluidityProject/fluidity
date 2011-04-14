/* Copyright (C) 2006 Imperial College London.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

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

#ifndef ERRORMEASURE_H
#define ERRORMEASURE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "vtk.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "MetricTensor.h"
#include "cinterfaces.h"

/** Used to form the metric tensor field used to direct mesh adaptivity.
 */
class ErrorMeasure{
 public:
  /// Constructor.
  ErrorMeasure();

  /// Destructor.
  ~ErrorMeasure();
  
  /// Add a field that contribures to the construction of the error metric tensor field.
  void add_field(std::string name, double error, bool relative, double sigma);
  
  /// Apply mesh gradation. Algorithm based on:
  /// Anisotropic mesh gradation control (2004), 
  /// Xiangrong Li, Jean-françois Remacle, Nicolas Chevaugeon,
  /// Mark S. Shephard in Thirteenth International Meshing Roundtable.
  void apply_gradation(double gradation);

  /// Add diagnostic fields to the vtkUnstructuredGrid
  void diagnostics();

  /// Get the expected number of elements after the mesh is adapted to the metric field.
  int get_expected_nelements();

  /// Get output.
  vtkUnstructuredGrid *get_output();

  /*
  /// Limit the aspect ratio of the desired edge lengths.
  int limit_aspect(double);
  */

  /// Set mesh information
  void set_input(vtkUnstructuredGrid *ug);

  /// Set the anisotropic maximum edge lengths allowed as
  /// a dim x dim symmetric tensor. If len==1 then is's the same
  /// limits everywhere. Otherwise this should be the same as the
  /// number of nodes in the mesh.
  void set_max_length(double *max_len, int len);

  /// Set homogenious maximum edge lengths.
  void set_max_length(double max_len);

  /// Set the maximum number of nodes permitted.
  void set_max_nodes(int max_nodes);

  /// Set the anisotropic minimum edge lengths allowed as
  /// a dim x dim symmetric tensor. If len==1 then is's the same
  /// limits everywhere. Otherwise this should be the same as the
  /// number of nodes in the mesh.
  void set_min_length(double *min_len, int len);

  /// Set homogenious minimum edge lengths.
  void set_min_length(double min_len);

  /// Turn off verbose mode
  void verbose_off();

  /// Turn on verbose mode
  void verbose_on();

 private:
  /// Abstraction of BLAS [ds]pev.
  int blas_spev(char jobz, char uplo, int N, const double ap[],
                double eigenvalues[], double eigenvectors[], int ldz,  double work[]);

  /// Abstraction of BLAS [ds]gemm.
  int blas_sgemm(char TRANSA, char TRANSB, int N, double A[], double B[], double C[]);

  /// Calculate Hessian of field.
  void get_hessian(std::string name);

  /// Gradient of field using SPR.
  int grad(std::string name);

  /// Cartesian distance between two points
  inline double length(size_t, size_t) const;
  
  /// Calculate volume of element
  // double calc_volume(int e) const;
 
  // Mesh
  vtkUnstructuredGrid *ug;
  size_t dim;

  std::vector< std::set<size_t> > NNList, NEList;
  std::set<size_t> boundary;

  bool verbose;

};

#endif
