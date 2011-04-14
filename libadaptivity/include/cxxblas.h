/* Copyright (C) 2006 Imperial College London and others.
 *
 *  Please see the AUTHORS file in the main source directory for a full list
 *  of copyright holders.
 *
 *  Dr Gerard J Gorman
 *  Applied Modelling and Computation Group
 *  Department of Earth Science and Engineering
 *  Imperial College London
 *
 *  g.gorman@imperial.ac.uk
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 *  USA
 */

#ifndef CXXBLAS_H
#define CXXBLAS_H

extern "C" {
  /* From BLAS: compute all the eigenvalues and, optionally,
 *      eigenvectors of a real symmetric matrix A in packed storage.
 *        */
  void sspev_(char *, char *, int *, const float [],  float [],  float [],  int *,  float [], int *);
  void dspev_(char *, char *, int *, const double [], double [], double [], int *, double [], int *);

  void sgemm_(const char *, const char *, const int *, const int *, const int *, 
              const float *, const float [], const int *, 
              const float [], const int *, const float *, float [], const int *);
  void dgemm_(const char *, const char *, const int *, const int *, const int *, 
              const double *, const double [], const int *, 
              const double [], const int *, const double *, double [], const int *);
}
#endif
