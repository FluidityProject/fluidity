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

#include "Fldecomp_Wrappers.h"

using namespace std;

namespace Fluidity{

  int ReadMesh(const string& filename, const string& meshType, vector<double>& x, int& dim,
               vector<int>& ENList, vector<int>& regionIds, int& nloc, vector<int>& SENList, vector<int>& boundaryIds, int& snloc){

    // Read mesh information so that memory can be allocated
    int filename_len;
    filename_len = filename.size();
    int meshType_len;
    meshType_len = meshType.size();
    int NElements, NNodes, SNElements;
    query_mesh(filename.c_str(), &filename_len, meshType.c_str(), &meshType_len,
                  &dim, &NNodes, &NElements, &nloc, &SNElements, &snloc);

    // Allocate memory to read the mesh
    double* x_handle = (double*)malloc(NNodes*dim*sizeof(double));
    int* ENList_handle = (int*)malloc(NElements*nloc*sizeof(int));
    int* regionIds_handle = (int*)malloc(NElements*sizeof(int));
    int* SENList_handle = (int*)malloc(SNElements*snloc*sizeof(int));
    int* boundaryIds_handle = (int*)malloc(SNElements*sizeof(int));
    if(x_handle == NULL or ENList_handle == NULL or regionIds_handle == NULL or SENList_handle == NULL or boundaryIds_handle == NULL){
      cerr << "Failed to allocate memory in ReadMesh" << endl;
      free(x_handle);
      free(ENList_handle);
      free(regionIds_handle);
      free(SENList_handle);
      free(boundaryIds_handle);
      return -1;
    }
    
    // Read the mesh
    int ret = read_mesh(filename.c_str(), &filename_len, meshType.c_str(), &meshType_len,
                        x_handle, &dim, &NNodes,
                        ENList_handle, regionIds_handle, &NElements, &nloc,
                        SENList_handle, boundaryIds_handle, &SNElements, &snloc);

    if(ret == 0){
      // Copy the mesh into STL data structures
      x.clear();  x.resize(NNodes * dim);
      for(int i = 0;i < NNodes * dim;i++){
        x[i] = x_handle[i];
      }
      ENList.clear();  ENList.resize(NElements * nloc);
      for(int i = 0;i < NElements * nloc;i++){
        ENList[i] = ENList_handle[i];
      }
      regionIds.clear();  regionIds.resize(NElements);
      for(int i = 0;i < NElements;i++){
        regionIds[i] = regionIds_handle[i];
      }
      SENList.clear();  SENList.resize(SNElements * snloc);
      for(int i = 0;i < SNElements * snloc;i++){
        SENList[i] = SENList_handle[i];
      }
      boundaryIds.clear();  boundaryIds.resize(SNElements);
      for(int i = 0;i < SNElements;i++){
        boundaryIds[i] = boundaryIds_handle[i];
      }
    }

    // Deallocate memory used to read the mesh
    free(x_handle);
    free(ENList_handle);
    free(regionIds_handle);
    free(SENList_handle);
    free(boundaryIds_handle);
    
    return ret;
  }

  int ReadMesh(const string& filename, const string& meshType, vector<double>& x, int& dim,
               vector<int>& ENList, vector<int>& regionIds, int& nloc, deque<vector<int> >& SENList, vector<int>& boundaryIds, int& snloc){
      
    vector<int> SENList_handle;
    int ret = ReadMesh(filename, meshType, x, dim,
                        ENList, regionIds, nloc, SENList_handle, boundaryIds, snloc);
      
    if(ret == 0){
    
      size_t SNElements = boundaryIds.size();
      assert(SENList_handle.size() == SNElements * snloc);
      SENList.clear();  SENList.resize(SNElements, vector<int>(snloc));
      for(size_t i = 0;i < SNElements;i++){
        for(int j = 0;j < snloc;j++){
          SENList[i][j] = SENList_handle[i * snloc + j];
        }
      }
    }
    
    return ret;
  }

  int WriteMesh(const string& filename, const string& meshType, const vector<double>& x, const int& dim,
                const vector<int>& ENList, const vector<int>& regionIds, const int& nloc, const vector<int>& SENList, const vector<int>& boundaryIds, const int& snloc){

    // Extract mesh information so that memory can be allocated
    int NNodes = x.size() / dim;
    int NElements = regionIds.size();
    assert((int)ENList.size() == NElements * nloc);
    int SNElements = boundaryIds.size();
    assert((int)SENList.size() == SNElements * snloc);
        
    // Allocate memory to write the mesh
    double* x_handle = (double*)malloc(NNodes*dim*sizeof(double));
    int* ENList_handle = (int*)malloc(NElements*nloc*sizeof(int));
    int* regionIds_handle = (int*)malloc(NElements*sizeof(int));
    int* SENList_handle = (int*)malloc(SNElements*snloc*sizeof(int));
    int* boundaryIds_handle = (int*)malloc(SNElements*sizeof(int));
    if(x_handle == NULL or ENList_handle == NULL or regionIds_handle == NULL or SENList_handle == NULL or boundaryIds_handle == NULL){
      cerr << "Failed to allocate memory in WriteMesh" << endl;
      free(x_handle);
      free(ENList_handle);
      free(regionIds_handle);
      free(SENList_handle);
      free(boundaryIds_handle);
      return -1;
    }

    // Copy the mesh out of STL data structures  
    for(int i = 0;i < NNodes * dim;i++){
      x_handle[i] = x[i];
    }
    for(int i = 0;i < NElements * nloc;i++){
      ENList_handle[i]=ENList[i];
    }
    for(int i = 0;i < NElements;i++){
      regionIds_handle[i] = regionIds[i];
    }
    for(int i = 0;i < SNElements * snloc;i++){
      SENList_handle[i] = SENList[i];
    }
    for(int i = 0;i < SNElements;i++){
      boundaryIds_handle[i] = boundaryIds[i];
    }

    // Write the mesh
    int filename_len;
    filename_len = filename.size();
    int meshType_len;
    meshType_len = meshType.size();
    int ret = write_mesh(filename.c_str(), &filename_len, meshType.c_str(), &meshType_len,
                         x_handle, &dim, &NNodes,
                         ENList_handle, regionIds_handle, &NElements, &nloc,
                         SENList_handle, boundaryIds_handle, &SNElements, &snloc);

    // Deallocate memory used to write the mesh
    free(x_handle);
    free(ENList_handle);
    free(regionIds_handle);
    free(SENList_handle);
    free(boundaryIds_handle);
    
    return ret;
  }

  int WriteMesh(const string& filename, const string& meshType, const vector<double>& x, const int& dim,
                const vector<int>& ENList, const vector<int>& regionIds, const int& nloc, const deque<vector<int> >& SENList, const vector<int>& boundaryIds, const int& snloc){

    size_t SNElements = boundaryIds.size();
    assert(SENList.size() == SNElements);
    for(size_t i = 0;i < SNElements;i++){
      assert((int)SENList[i].size() == snloc);
    }
        
    vector<int> SENList_handle(SENList.size() * snloc);
    for(size_t i = 0;i < SNElements;i++){
      for(int j = 0;j < snloc;j++){
        SENList_handle[i * snloc + j] = SENList[i][j];
      }
    }
    
    return WriteMesh(filename, meshType, x, dim,
                      ENList, regionIds, nloc, SENList_handle, boundaryIds, snloc);
  }
  
}
