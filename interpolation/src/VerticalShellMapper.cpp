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
#include "confdefs.h"

#include <string.h>

#include "VerticalShellMapper.h"

using namespace std;
using namespace SpatialIndex;

class MyVisitor : public IVisitor{
public:
  MyVisitor(){
    nearest.first = -1;
    for(size_t i=0;i<3;i++)
      nearest_shape[i] = -1.0;
  }
  
  int GetHit(double *shape) const{
    for(size_t i=0;i<3;i++)
      shape[i] = nearest_shape[i];
    
    return nearest.first;
  }

  void SetCoordinates(const double *_sx, const double *_sy){
    sx = _sx;
    sy = _sy;
  }

  void SetPoint(double x, double y){
    find_x = x;
    find_y = y;
  }

  void visitNode(const INode& n) {}
  
  void visitData(const IData& d){
    int id=d.getIdentifier();
    
    int *nodes;
#ifdef EXPERIMENTAL_SPATIALINDEX
    uint32_t len;
#else
    size_t len;
#endif
    d.getData(len, reinterpret_cast<byte**>(&nodes));

    double shape[3];
    double io = inside_element(find_x, find_y, id, nodes[0], nodes[1], nodes[2], shape);
    
    if(nearest.first==-1){
      nearest.first = id;
      nearest.second = io;
      for(size_t i=0;i<3;i++)
        nearest_shape[i] = shape[i];
    }else{
      if(io>nearest.second){
        nearest.first = id;
        nearest.second = io;
        for(size_t i=0;i<3;i++)
          nearest_shape[i] = shape[i];
      }
    }
    
    delete [] nodes;
  }
  
  void visitData(std::vector<const IData*>& v) {}
  
private:
  double inside_element(double x, double y, int id, int n0, int n1, int n2, double *shape){
    if(normals.find(id)==normals.end())
      normals[id] = normal(sx[n0], sy[n0], sx[n1], sy[n1], sx[n2], sy[n2]);
    
    double Normal = normals[id];
    
    double N=normal(sx[n0], sy[n0],
                    sx[n1], sy[n1],
                    x, y);
    shape[2] = N/Normal;
    
    N = normal(sx[n1], sy[n1],
               sx[n2], sy[n2],
               x, y);
    shape[0] = N/Normal;
    
    N = normal(sx[n2], sy[n2],
               sx[n0], sy[n0],
               x, y);
    shape[1] = N/Normal;
    
    return min(shape[0], min(shape[1], shape[2]));
  }
  
  double normal(double x0, double y0,
                double x1, double y1,
                double x2, double y2) const{
    double AB[2];
    AB[0] = x1 - x0;
    AB[1] = y1 - y0;
    
    double AC[2];
    AC[0] = x2 - x0;
    AC[1] = y2 - y0;
    
    // normal = AB X AC
    return (AB[0]*AC[1] - AC[0]*AB[1]);
  }
  
  const double *sx, *sy;
  map<int, double> normals;
  pair<int, double> nearest;
  double nearest_shape[3];
  double find_x, find_y;
};

VerticalShellMapper::VerticalShellMapper(){
  verbose = false;

  spherical_geometry = false;

  data_block_size = 3*sizeof(int);

  // Create storage managers
  storage_manager_north = StorageManager::createNewMemoryStorageManager();
  storage_manager_south = StorageManager::createNewMemoryStorageManager();
  storage_manager = StorageManager::createNewMemoryStorageManager();
  
  // R-Tree parameters
  double fillFactor = 0.7;
  unsigned long indexCapacity = 100;
  unsigned long leafCapacity = 100;
  unsigned long dimension = 2;
  RTree::RTreeVariant variant = RTree::RV_RSTAR;
  
  // create R-trees
  indexId_north=1;
  rtree_north = RTree::createNewRTree(*storage_manager_north, fillFactor, indexCapacity,
                                      leafCapacity, dimension, variant, indexId_north);

  indexId_south=2;
  rtree_south = RTree::createNewRTree(*storage_manager_south, fillFactor, indexCapacity,
                                      leafCapacity, dimension, variant, indexId_south);

  indexId=3;
  rtree = RTree::createNewRTree(*storage_manager, fillFactor, indexCapacity,
                                leafCapacity, dimension, variant, indexId);
  // rtree = RTree::loadRTree(*storage_manager, 1);
}

VerticalShellMapper::~VerticalShellMapper(){
  delete rtree_north;
  delete rtree_south;
  delete rtree;

  delete storage_manager_north;
  delete storage_manager_south;
  delete storage_manager;
}

void VerticalShellMapper::DumpVTU(const char *basename) const{
  if(verbose)
    cout<<"void VerticalShellMapper::DumpVTU("<<basename<<") const\n";
  
#ifdef HAVE_VTK
  if(Triangles_north.size()>0){ // north
    string filename=string(basename)+"_north.vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();
    
    // Point definitions
    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    const size_t ncnt = sx.size();
    for(size_t i=0; i<ncnt; i++){
      newPts->InsertNextPoint(sx[i], sy[i], 0.0);
    }
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();
    
    vtkIdType Cell[3];
    const size_t ecnt = Triangle_ids_north.size();
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<3; j++){
        assert(Triangles_north[i*3+j]>=0);
        assert(Triangles_north[i*3+j]<(int)ncnt);
        Cell[j] = Triangles_north[i*3+j];
      }
      ug->InsertNextCell(VTK_TRIANGLE, 3, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }

  if(Triangles_south.size()>0){ // south
    string filename=string(basename)+"_south.vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();
    
    // Point definitions
    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    const size_t ncnt = sx.size();
    for(size_t i=0; i<ncnt; i++)
      newPts->InsertNextPoint(sx[i], sy[i], 0.0);
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();
    
    vtkIdType Cell[3];
    const size_t ecnt = Triangle_ids_south.size();
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<3; j++){
        assert(Triangles_south[i*3+j]>=0);
        assert(Triangles_south[i*3+j]<(int)ecnt);
        Cell[j] = Triangles_south[i*3+j];
      }
      ug->InsertNextCell(VTK_TRIANGLE, 3, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }

  if(Triangles.size()>0){ // plane
    string filename=string(basename)+"_plane.vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();
    
    // Point definitions
    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    const size_t ncnt = sx.size();
    for(size_t i=0; i<ncnt; i++)
      newPts->InsertNextPoint(sx[i], sy[i], 0.0);
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();
    
    vtkIdType Cell[3];
    const size_t ecnt = Triangles.size()/3;
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<3; j++){
        assert(Triangles[i*3+j]>=0);
        assert(Triangles[i*3+j]<(int)ecnt);
        Cell[j] = Triangles[i* 3+j];
      }
      ug->InsertNextCell(VTK_TRIANGLE, 3, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }
#endif
}

// Return value is the element number which contains the coordinate
// under consideration. -1 is returned if no element is found. This
// should only happen if the element is actually on another processes.
int VerticalShellMapper::Find(double x, double y, double z, double *shape){
  int tri_id=-1;
  double pos[2];

  if(spherical_geometry){

    bool north;
    if(z>=0){
      projection::cartesian2stereographic(0.0, 90.0, x, y, z, pos[0], pos[1]);
      north = true;
      if(verbose)
        cout<<"Searching northern hemisphere\n";
    }else{
      projection::cartesian2stereographic(0.0,-90.0, x, y, z, pos[0], pos[1]);
      north = false;
      if(verbose)
        cout<<"Searching southern hemisphere\n";
    }
    
    Point p = Point(pos, 2);
    
    MyVisitor v;
    v.SetCoordinates(&(sx[0]), &(sy[0]));
    v.SetPoint(pos[0], pos[1]);
    
    if(north){
      rtree_north->pointLocationQuery(p, v);
    }else{
      rtree_south->pointLocationQuery(p, v);
    }

    tri_id = v.GetHit(shape);

  }else{
    pos[0] = x;
    pos[1] = y;
    
    Point p = Point(pos, 2);
    
    MyVisitor v;
    v.SetCoordinates(&(sx[0]), &(sy[0]));
    v.SetPoint(pos[0], pos[1]);
    
    rtree->pointLocationQuery(p, v);

    tri_id = v.GetHit(shape);
  }
  
  return tri_id;
}

bool VerticalShellMapper::IsParallel() const{
  if(verbose)
    cout<<"bool VerticalShellMapper::IsParallel() const\n";
  
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    return (MPI::COMM_WORLD.Get_size()>1);
  }else{
    return false;
  }
#else
  return false;
#endif
}

// Return the process rank which has target element and the id of the
// said element. -1 is returned if no element is found and this should
// be treated as an error - possibly the free surface blown up.
int VerticalShellMapper::pFind(const double *x, const double *y, const double *z, int nnodes,
                               int *tid, double *shape, int *host){
  size_t myrank = 0;
  size_t NProcs = 1;
  
#ifdef HAVE_MPI
  if(IsParallel()){
    myrank = MPI::COMM_WORLD.Get_rank();
    NProcs = MPI::COMM_WORLD.Get_size();
  }
#endif
  
  vector<double> xyz;
  vector<int> lut;
  for(size_t i=0;i<(size_t)nnodes;i++){
    if(z==NULL)
      tid[i] = Find(x[i], y[i], 0.0, shape+i*3);
    else
      tid[i] = Find(x[i], y[i], z[i], shape+i*3);
    
    double min_weight = min(shape[0], min(shape[1], shape[2]));
    if(min_weight<0.0){ // .. continue on from here
      xyz.push_back(x[i]);
      xyz.push_back(y[i]);
      if(z==NULL)
        xyz.push_back(0.0);
      else
        xyz.push_back(z[i]);
      
      lut.push_back(i);
    }else{
      host[i] = myrank;
    }
  }

#ifdef HAVE_MPI
  if(IsParallel()){
    int missing=xyz.size();
    int total_missing=0;
    MPI::COMM_WORLD.Allreduce(&missing, &total_missing, 1, MPI::INT, MPI::SUM);
    if(total_missing>0){
      for(size_t i=0;i<NProcs;i++){
        MPI::COMM_WORLD.Bcast(&missing, 1, MPI::INT, i);
        if(myrank!=i){
          xyz.resize(missing);
        }
        MPI::COMM_WORLD.Bcast(&(xyz[0]), missing, MPI::DOUBLE, i);
        
        if(myrank!=i){
          missing /= 3;
          vector<int> found_tid;
          vector<double> missing_shape(3), found_shape;
          for(size_t j=0;j<(size_t)missing;j++){
            int missing_tid = Find(xyz[j*3], xyz[j*3+1], xyz[j*3+2], &(missing_shape[0]));
            if(missing_tid>=0){
              found_tid.push_back(j);
              found_tid.push_back(missing_tid);
              found_shape.insert(found_shape.end(), missing_shape.begin(), missing_shape.end());
            }
          }
          
          if(found_tid.size()>0){
            MPI::COMM_WORLD.Send(&(found_tid[0]), found_tid.size(), MPI::INT, i, 33);
            MPI::COMM_WORLD.Send(&(found_shape[0]), found_shape.size(), MPI::DOUBLE, i, 34);
          }else{
            MPI::COMM_WORLD.Send(NULL, 0, MPI::INT, i, 33);
          }
        }else{
          for(size_t j=0;j<NProcs;j++){
            if(i==j)
              continue;
            
            MPI::Status recvStatus;
            MPI::COMM_WORLD.Probe(j, 33, recvStatus);  
            int mesg_size = recvStatus.Get_count(MPI::INT);
            int found_cnt = mesg_size/2;
            
            vector<int> buffer_int(found_cnt*2);
            MPI::COMM_WORLD.Recv(&(buffer_int[0]), found_cnt*2, MPI::INT, j, 33);
            
            vector<double> buffer_real(3*found_cnt);
            if(found_cnt)
              MPI::COMM_WORLD.Recv(&(buffer_real[0]), 3*found_cnt, MPI::DOUBLE, j, 34);
            
            for(int k=0;k<found_cnt;k++){
              tid[lut[buffer_int[k*2]]] = buffer_int[k*2+1];
              for(int l=0;l<3;l++)
                shape[lut[k*2]+l] = buffer_real[k*3+l];
              host[lut[buffer_int[k*2]]] = j;
            }
          }
        }
      }
    }
  }
#endif
  
  for(size_t i=0;i<(size_t)nnodes;i++)
    if(tid[i]<0)
      return -1;

  return 0;
}


void VerticalShellMapper::SetShell(const double *x, const double *y, const double *z,
                                   const int *SurfaceTriangles, int NTriangles, int zero_index){
  if(spherical_geometry){
    // Create two sets of elements - those that has nodes in the
    // northern helisphere and those with elements in the southern
    // hemisphere. Note that elements may appear in both sets and their
    // nodes will be duplicated. Each hemisphere will be mapped with a
    // stereographic projection using the opposing pole as the
    // projection point which goes to infinity.
    
    set<int> elements_north, elements_south;
    map<int, int> renumber_north, renumber_south;
    int mxnid=0;
    for(size_t i=0;i<(size_t)NTriangles;i++){
      bool north=false, south=false;
      
      for(size_t j=0;j<3;j++){
        int nid = SurfaceTriangles[i*3+j]-zero_index;
        if(z[nid]>=0)
          north = true;
        else
          south = true;
      }
      
      if(north){
        elements_north.insert(i);
        for(size_t j=0;j<3;j++){
          int nid = SurfaceTriangles[i*3+j]-zero_index;
          renumber_north[nid] = -1;
          mxnid=max(nid, mxnid);
        }
      }
      
      if(south){    
        elements_south.insert(i);
        for(size_t j=0;j<3;j++){
          int nid = SurfaceTriangles[i*3+j]-zero_index;
          renumber_south[nid] = -1;
        }
      }
    }
    
    for(int i=0;i<mxnid;i++)
      renumber_north[i] = -1;

    // Finalise renumbering
    int pos=0;
    for(map<int, int>::iterator it=renumber_north.begin(); it!=renumber_north.end(); ++it){
      it->second = pos++;
    }
    
    for(map<int, int>::iterator it=renumber_south.begin(); it!=renumber_south.end(); ++it){
      it->second = pos++;
    }
    
    // Store renumbered elements
    for(set<int>::const_iterator it=elements_north.begin(); it!=elements_north.end(); ++it){
      Triangle_ids_north.push_back(*it);
      for(size_t j=0;j<3;j++){
        int nid = SurfaceTriangles[(*it)*3+j]-zero_index;
        Triangles_north.push_back(renumber_north[nid]);
      }
    }
    
    for(set<int>::const_iterator it=elements_south.begin(); it!=elements_south.end(); ++it){
      Triangle_ids_south.push_back(*it);
      for(size_t j=0;j<3;j++){
        int nid = SurfaceTriangles[(*it)*3+j]-zero_index;
        Triangles_south.push_back(renumber_south[nid]);
      }
    }
    
    // Project and store all the valid points
    for(map<int, int>::const_iterator it=renumber_north.begin(); it!=renumber_north.end(); ++it){
      double xx, yy;
      projection::cartesian2stereographic(0.0, 90.0, x[it->first], y[it->first], z[it->first], xx, yy);
      sx.push_back(xx);
      sy.push_back(yy);
    }
    
    for(map<int, int>::const_iterator it=renumber_south.begin(); it!=renumber_south.end(); ++it){
      double xx, yy;
      projection::cartesian2stereographic(0.0, -90.0, x[it->first], y[it->first], z[it->first], xx, yy);
      sx.push_back(xx);
      sy.push_back(yy);
    }
    
    // Fill north RTree
    double plow[2], phigh[2];
    for(size_t i=0;i<Triangle_ids_north.size();i++){
      plow[0] = sx[Triangles_north[i*3]]; phigh[0] = plow[0];
      plow[1] = sy[Triangles_north[i*3]]; phigh[1] = plow[1];
      for(size_t j=1;j<3;j++){
        plow[0]  = min(plow[0],  sx[Triangles_north[i*3+j]]);
        phigh[0] = max(phigh[0], sx[Triangles_north[i*3+j]]);
        
        plow[1]  = min(plow[1],  sy[Triangles_north[i*3+j]]);
        phigh[1] = max(phigh[1], sy[Triangles_north[i*3+j]]);
      }
      
      double dx = (phigh[0] - plow[0])*0.01;
      double dy = (phigh[1] - plow[1])*0.01;

      plow[0] = plow[0] - dx;
      phigh[0] = phigh[0] + dx;

      plow[1] = plow[1] - dy;
      phigh[1] = phigh[1] + dy;
        
      Region r = Region(plow, phigh, 2);
      
      rtree_north->insertData(data_block_size, reinterpret_cast<const byte*>(&(Triangles_north[i*3])),
                              r, Triangle_ids_north[i]);
    }
    
    // Fill south RTree
    for(size_t i=0;i<Triangle_ids_south.size();i++){
      plow[0] = sx[Triangles_south[i*3]]; phigh[0] = plow[0];
      plow[1] = sy[Triangles_south[i*3]]; phigh[1] = plow[1];
      for(size_t j=1;j<3;j++){
        plow[0]  = min(plow[0],  sx[Triangles_south[i*3+j]]);
        phigh[0] = max(phigh[0], sx[Triangles_south[i*3+j]]);
        
        plow[1]  = min(plow[1],  sy[Triangles_south[i*3+j]]);
        phigh[1] = max(phigh[1], sy[Triangles_south[i*3+j]]);
      }
      
      double dx = (phigh[0] - plow[0])*0.01;
      double dy = (phigh[1] - plow[1])*0.01;
      
      plow[0] = plow[0] - dx;
      phigh[0] = phigh[0] + dx;
      
      plow[1] = plow[1] - dy;
      phigh[1] = phigh[1] + dy;
      
      Region r = Region(plow, phigh, 2);
      
      rtree_south->insertData(data_block_size, reinterpret_cast<const byte*>(&(Triangles_south[i*3])),
                              r, Triangle_ids_south[i]);
    }
  }else{
    // Store renumbered elements
    map<int, int> renumber;
    for(size_t i=0;i<(size_t)NTriangles;i++){
      for(size_t j=0;j<3;j++){
        int nid = SurfaceTriangles[i*3+j]-zero_index;
        renumber[nid] = -1;
        Triangles.push_back(nid);
      }
    }
    
    int pos=0;
    for(map<int, int>::iterator it=renumber.begin(); it!=renumber.end(); ++it){
      it->second = pos++;
    }
    
    for(vector<int>::iterator it=Triangles.begin(); it!=Triangles.end(); ++it)
      *it = renumber[*it];
  
    // Store all the valid points
    for(map<int, int>::const_iterator it=renumber.begin(); it!=renumber.end(); ++it){
      sx.push_back(x[it->first]);
      sy.push_back(y[it->first]);
    }
    
    // Fill RTree
    double plow[2], phigh[2];
    for(size_t i=0;i<(size_t)NTriangles;i++){
      int nid = Triangles[i*3];
      
      plow[0] = sx[nid]; phigh[0] = plow[0];
      plow[1] = sy[nid]; phigh[1] = plow[1];
      for(size_t j=1;j<3;j++){
        int nid = Triangles[i*3+j];
        
        plow[0]  = min(plow[0],  sx[nid]);
        phigh[0] = max(phigh[0], sx[nid]);
        
        plow[1]  = min(plow[1],  sy[nid]);
        phigh[1] = max(phigh[1], sy[nid]);
      }
  
      // Add tolerance
      double dx=0.01*(phigh[0] - plow[0]);
      phigh[0]+=dx; plow[0]-=dx;

      double dy=0.01*(phigh[1] - plow[1]);
      phigh[1]+=dy; plow[1]-=dy;

      Region r = Region(plow, phigh, 2);
      
      rtree->insertData(data_block_size, reinterpret_cast<const byte*>(&(Triangles[i*3])), r, i);
    }
  }  
}

void VerticalShellMapper::SetSphericalGeometry(bool stat){
  spherical_geometry = stat;
}

extern "C"{
#define verticalshellmapper_find_fc F77_FUNC_(verticalshellmapper_find, VERTICALSHELLMAPPER_FIND)
  void verticalshellmapper_find_fc(const double x[], const double y[], const double z[], const int *nnodes,
                                   const int senlist[], const int *ntri, const int *fortran_zero_index,
                                   const int *flat_earth, int tri_ids[], double shape_fxn[]){
    if((*ntri)<1){
      cerr<<"ERROR: verticalshellmapper_find() is called but "
          <<"no surface elements are available. This suggests "
          <<"that the top elements have not been indicated.\n";
    }
    
    VerticalShellMapper locator;
    locator.SetSphericalGeometry((*flat_earth)==0);
    locator.SetShell(x, y, z, senlist, *ntri, *fortran_zero_index);
    
    
    vector<int> host(*nnodes);
    locator.pFind(&(x[0]), &(y[0]), &(z[0]), *nnodes,
                  tri_ids, shape_fxn, &(host[0]));
    
    if(*fortran_zero_index)
      for(size_t i=0;i<(size_t)(*nnodes);i++)
        tri_ids[i] = tri_ids[i] + (*fortran_zero_index);
    
    return;
  }
  
#define mesh2npoints_2d_find_fc F77_FUNC_(mesh2npoints_2d_find, MESH2NPOINTS_2D_FIND)
  void mesh2npoints_2d_find_fc(const double x_src[], const double y_src[], const int *nnodes_src,
                               const int senlist_src[], const int *ntri_src, const int *fortran_zero_index,
                               const double x_dest[], const double y_dest[], const int *nnodes_dest,
                               int tri_ids[], double shape_fxn[]){
    if((*ntri_src)<1){
      cerr<<"ERROR:  mesh2npoints_2d_find() is called but "
          <<"no 2d elements.\n";
    }
    
    VerticalShellMapper locator;
    locator.SetSphericalGeometry(false);
    locator.SetShell(x_src, y_src, NULL, senlist_src, *ntri_src, *fortran_zero_index);

    vector<int> host(*nnodes_dest);
    locator.pFind(&(x_dest[0]), &(y_dest[0]), NULL, *nnodes_dest,
                  tri_ids, shape_fxn, &(host[0]));

    if(*fortran_zero_index)
      for(size_t i=0;i<(size_t)(*nnodes_dest);i++)
        tri_ids[i] = tri_ids[i] + (*fortran_zero_index);

  }
}

#ifdef TEST_VERTICALSHELLMAPPER

int main(){
#ifdef HAVE_VTK
  char filename[]="globe_surface.vtu";
  vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(filename);
  reader->Update();
  vtkUnstructuredGrid *ug=reader->GetOutput();
  
  vtkIdType ncnt = ug->GetNumberOfPoints();
  vector<double> x(ncnt), y(ncnt), z(ncnt);
  for(vtkIdType i=0;i<ncnt;i++){
    double *r = ug->GetPoint(i);
    x[i] = r[0];
    y[i] = r[1];
    z[i] = r[2];
  }
  
  vtkIdType ecnt = ug->GetNumberOfCells();
  vector<int> ENList(ecnt*3);
  for(vtkIdType i=0;i<ecnt;i++){
    for(size_t j=0;j<3;j++){
      ENList[i*3+j] = ug->GetCell(i)->GetPointId(j);
    }
  }
  
  VerticalShellMapper vertical;
  vertical.SetSphericalGeometry(true);
  vertical.SetShell(&(x[0]), &(y[0]), &(z[0]), &(ENList[0]), ecnt, 0);  
  vertical.DumpVTU("projection");

  double shape[3];
  for(vtkIdType i=0;i<ncnt;i++){
    cout<<"Looking for element containing node: "<<i<<", ("<<x[i]<<", "<<y[i]<<", "<<z[i]<<")"<<endl;
    int eid = vertical.Find(x[i], y[i], z[i], shape);
    cout<<"Found element: "<<eid<<". Contains:\n";
    for(size_t j=0;j<3;j++){
      int nid = ug->GetCell(eid)->GetPointId(j);
      cout<<nid<<":\t"<<x[nid]<<" "<<y[nid]<<" "<<z[nid]<<endl;
    }
  }
  
#else
  cerr<<"ERROR: no vtk support compiled\n";
#endif
  
  return 0;
}
#endif
