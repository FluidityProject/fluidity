/* Copyright (C) 2006 Imperial College London and others.

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

#ifndef METRICTENSOR_H
#define METRICTENSOR_H

#include "confdefs.h"

#include "vtk.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <float.h>

#include "cxxblas.h"

#ifdef FLT_SMALL
#error FLT_SMALL is already defined
#endif
#define FLT_SMALL (1.0E-12)

/** Performs basic operations on metric tensors.
 */
class MetricTensor{
 public:
  /// Constructor.
  MetricTensor(int dimension);

  /// Constructor.
  MetricTensor(double m11,
               double m12, double m22);

  /// Constructor.
  MetricTensor(double m11,
               double m12, double m22,
               double m13, double m23, double m33);

  /// Constructor - T stores the bottom triangle of the tensor.
  MetricTensor(int dim, const double *T);

  /// Copy constructor.
  MetricTensor(const MetricTensor&);

  /// Overloaded assignment operator.
  const MetricTensor& operator=(const MetricTensor &); 

  /// Overloaded write operator.
  friend std::ostream &operator<< (std::ostream&, const MetricTensor&);

  /// Aspect ratio between the maximum and minimum edge lengths
  /// required by the metric.
  double aspect_ratio() const;

  /// Average desired edge length
  double average_length() const;
  
  /// Calculates metric tensor that circumscribes this tensor and
  /// input, M. This limits the minimum edge lengths anisotropically.
  int circumscribe(const MetricTensor& M);

  /// Determinant of tensor
  double det() const;

  /// Calculate eivenvalue decomposition of metric tensor.
  /// @param D Eigenvalues
  /// @param V Eigenvectors.
  int eigen_decomp(double *D, double *V) const;

  /// Construct metric tensor from eigenvalue decomposition.
  /// @param D Eigenvalues
  /// @param V Eigenvectors.
  int eigen_undecomp(const double *D, const double *V);

  /// Calculate the inner product between this tensor and M2
  /// @returns Inner product.
  MetricTensor dot(const MetricTensor &M2);

  /// Pass back a copy of the bottom triangle of metric tensor into M (must be allocated by caller).
  void get_metric(double *M) const;

  /// Pass back a copy of the full metric tensor into M (must be allocated by caller).
  void get_metric2(double *M) const;

  /// Calculates metric tensor that inscribes this tensor and
  /// input, M. This limits the maximum edge lengths anisotropically.
  int inscribe(const MetricTensor& M);

  /// Limit the maximum edge length.
  int limit_max_size(const double *max_d);

  /// Limit the minimum edge length.
  int limit_min_size(const double *min_d);

  /// Calculate inverse of metric tensor.
  MetricTensor inv() const;

  /// True if metric is null.
  bool is_null() const;

  /// Limit the aspect ratio of the desired edge lengths.
  int limit_aspect(double);

  /// Scale the metric tensor.
  /// @param s Scale factor
  void scale(double s);

  /// Superimposes two metrics.
  int superimpose(const MetricTensor& M);

  /// Disable verbose mode.
  void verbose_off();
  
  /// Enable verbose mode.
  void verbose_on();
  
  /// Dump a copy of the metric as vtkPolyData. Must be configured with VTK support.
  int write_vtk(const char *filename) const;

 private:
  size_t dim, t_size;
  std::vector<double> metric;
  static bool verbose;

  /// BLAS [DS]PEV abstraction
  int blas_spev(char jobz, char uplo, int N, const double ap[],  
		double eigenvalues[],  double eigenvectors[],  int ldz,  double work[]) const;

  /// Cofactor of index i, j
  double cofactor(int i, int j) const;

  /// Calculate eivenvalue decomposition of tensor T.
  /// @param T tensor
  /// @param D Eigenvalues
  /// @param V Eigenvectors.
  int eigen_decomp(const double *T, double *D, double *V) const;

  /// Construct metric tensor from eigenvalue decomposition.
  /// @param T Tensor
  /// @param D Eigenvalues
  /// @param V Eigenvectors.
  int eigen_undecomp(const double *D, const double *V, double *T) const;


  /// In-place inversion of tensor.
  int inv(double&, double&,
	  double&, double&) const;

  /// In-place inversion of tensor.
  int inv(double&, double&, double&,
          double&, double&, double&,
	  double&, double&, double&) const;
  
  /// Convert i, j index into lower-triangle index
  int lookup(size_t i, size_t j) const;

  /// Map self to the same space where M is the unit sphere
  int map(const MetricTensor& M);

  int MatrixDotMatrix(double *A, bool aT, double *B, bool bT, double *C) const;
  void print_eigen(double *, double*) const;

  /// Unmap M back from unit metric space.
  int unmap(const MetricTensor& M);
};

#endif
