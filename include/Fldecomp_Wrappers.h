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

#ifndef FLDECOMP_WRAPPERS_H
#define FLDECOMP_WRAPPERS_H

#include <cassert>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"

namespace Fluidity{

  int ReadMesh(const std::string& filename, const std::string& meshType, std::vector<double>& x, int& dim,
                std::vector<int>& ENList, std::vector<int>& regionIds, int& nloc, std::vector<int>& SENList, std::vector<int>& boundaryIds, int& snloc);

  int ReadMesh(const std::string& filename, const std::string& meshType, std::vector<double>& x, int& dim,
                std::vector<int>& ENList, std::vector<int>& regionIds, int& nloc, std::deque<std::vector<int> >& SENList, std::vector<int>& boundaryIds, int& snloc);

  int WriteMesh(const std::string& filename, const std::string& meshType, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                 const std::vector<int>& ENList, const std::vector<int>& regionIds, const int& nloc, const std::vector<int>& SENList, const std::vector<int>& boundaryIds, const int& snloc);
                 
  int WriteMesh(const std::string& filename, const std::string& meshType, const std::vector<double>& x, const int& dim,
                 const std::vector<int>& ENList, const std::vector<int>& regionIds, const int& nloc, const std::deque<std::vector<int> >& SENList, const std::vector<int>& boundaryIds, const int& snloc);
}

extern "C"{

#define query_mesh F77_FUNC(query_mesh, QUERY_MESH)
  void query_mesh(const char* filename, const int* filename_len, const char* meshtype, const int* meshtype_len,
                  int* dim, int* nnodes, int* nelements, int* nloc, int* snelements, int* snloc);

#define read_mesh F77_FUNC(read_mesh, READ_MESH)
  int read_mesh(const char* filename, const int* filename_len, const char* meshtype, const int* meshtype_len,
                double* x, const int* dim, const int* nnodes,
                int* enlist, int* region_ids, const int* nelements, const int* nloc,
                int* senlist, int* boundary_ids, const int* snelements, const int* snloc);

#define write_mesh F77_FUNC(write_mesh, WRITE_MESH)
  int write_mesh(const char* filename, const int* filename_len, const char* meshtype, const int* meshtype_len,
                 const double* x, const int* dim, const int* nnodes,
                 const int* enlist, const int* region_ids, const int* nelements, const int* nloc,
                 const int* senlist, const int* boundary_ids, const int* snelements, const int* snloc);

}
                 
#endif
