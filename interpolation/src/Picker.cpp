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

#include "Picker.h"

using namespace std;

PickerMeshDataStream::PickerMeshDataStream(const double *x, const double *y, const double *z,
                                           const int *enlist, int zero_index, int nelements, int loc) : m_pNext(0){
  init(x, y, z, enlist, zero_index, nelements, loc, NULL);
}

PickerMeshDataStream::PickerMeshDataStream(const double *x, const double *y, const double *z,
                                           const int *enlist, int zero_index, int nelements, int loc,
                                           const int *live_set) : m_pNext(0){
  init(x, y, z, enlist, zero_index, nelements, loc, live_set);
}

PickerMeshDataStream::PickerMeshDataStream(const double *x, const double *y,
                                           const int *enlist, int zero_index, int nelements, int loc) : m_pNext(0){
  init(x, y, NULL, enlist, zero_index, nelements, loc, NULL);
}

PickerMeshDataStream::PickerMeshDataStream(const double *x, const double *y,
                                           const int *enlist, int zero_index, int nelements, int loc,
                                           const int *live_set) : m_pNext(0){
  init(x, y, NULL, enlist, zero_index, nelements, loc, live_set);
}

void PickerMeshDataStream::init(const double *x, const double *y, const double *z,
                                const int *enlist, int zero_index, int nelements, int loc,
                                const int *live_set){

  this->zero_index = zero_index;

  assert(x);
  assert(y);
  if(z)
    dim = 3;
  else
    dim = 2;

  assert(enlist);
  assert(nelements >= 0);
#ifdef DDEBUG
  if(nelements > 0){
    assert(loc >= 1);
  }
#endif

  this->x = x;
  this->y = y;
  this->z = z;
  
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  
  this->enlist = enlist;
  this->nelements = nelements;
  this->loc = loc;
  this->live_set = live_set;

  data_block_size = loc*sizeof(int);

  index = 0;

  readNextEntry();
}

PickerMeshDataStream::~PickerMeshDataStream(){
  if (m_pNext != 0) delete m_pNext;
}

SpatialIndex::IData* PickerMeshDataStream::getNext(){
  if (m_pNext == 0) return 0;
  
  SpatialIndex::RTree::Data* ret = m_pNext;
  m_pNext = 0;
  readNextEntry();
  return ret;
}

void PickerMeshDataStream::readNextEntry(){
  if(index >= nelements){
    return;
  }

  int eid;
  if(live_set)
    eid = live_set[index];
  else
    eid = index;
  
  index++;

  vector<int> node(loc);
  double high[3], low[3];

  node[0] = enlist[loc*eid] - zero_index;
  for(int j=0;j<dim;j++){
    low[j]  = pos[j][node[0]];
    high[j] = pos[j][node[0]];
  }
  
  for(int i=1;i<loc;i++){
    node[i] = enlist[loc*eid+i] - zero_index;
    for(int j=0;j<dim;j++){
      low[j]  = min(low[j],  pos[j][node[i]]);
      high[j] = max(high[j], pos[j][node[i]]);
    }
  }
  
  // Expand the bounding box, to help prevent missed queries
  for(int i=0;i<dim;i++){
    double dx = high[i] - low[i];
    low[i]  -= 0.01*dx;
    high[i] += 0.01*dx;
  }
    
  SpatialIndex::Region r(low, high, dim);
  m_pNext = new SpatialIndex::RTree::Data(data_block_size, reinterpret_cast<byte*>(&(node[0])), r, eid);

  return;
}

bool PickerMeshDataStream::hasNext(){
  return (m_pNext != 0);
}

#ifdef EXPERIMENTAL_SPATIALINDEX
uint32_t PickerMeshDataStream::size(){
#else
size_t PickerMeshDataStream::size(){
#endif
  return nelements;
}

void PickerMeshDataStream::rewind(){
  if (m_pNext != 0){
    delete m_pNext;
    m_pNext = 0;
  }
  
  index = 0;
  readNextEntry();
  return;
}

class SIVisitor : public SpatialIndex::IVisitor{
public:
  SIVisitor(unsigned long _element_type, unsigned long _dimension):
    TRI(5), QUAD(9), TETRA(10), HEX(12){
    
    element_type = _element_type;
    if(element_type==TRI)
      nloc = 3;
    else if(element_type==QUAD)
      nloc = 4;
    else if(element_type==TETRA)
      nloc = 4;
    else if(element_type==HEX)
      nloc = 8;
    else{
      cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): unknown element type in SIVisitor.\n";
      exit(-1);
    }
    nearest_shape.resize(nloc);
    for(size_t i=0;i<nloc;i++)
      nearest_shape[i] = -1.0;
    nearest.first = -1;
    
    dimension = _dimension;
    if(((element_type==TRI)||(element_type==QUAD))&&(dimension==3)){
      cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): unsupported element_type("<<element_type
          <<"), dimension("<<dimension<<") pair.\n";
      exit(-1);
    }
  }
  
  int GetHit(double *shape) const{
    for(size_t i=0;i<nloc;i++){
      shape[i] = nearest_shape[i];
    }
    
    return nearest.first;
  }

  void SetCoordinates(const double *_sx, const double *_sy){
    sx = _sx;
    sy = _sy;
  }

  void SetCoordinates(const double *_sx, const double *_sy, const double *_sz){
    sx = _sx;
    sy = _sy;
    sz = _sz;
  }

  void SetPoint(const double* point){
    find_x = point[0];
    find_y = point[1];
    if(dimension>2)
      find_z = point[2];
  }

  void visitNode(const SpatialIndex::INode& n) {}
  
  void visitData(const SpatialIndex::IData& d){
    int id=d.getIdentifier();
    
    int *nodes;
#ifdef EXPERIMENTAL_SPATIALINDEX
    uint32_t len;
#else
    size_t len;
#endif
    d.getData(len, reinterpret_cast<byte**>(&nodes));

    vector<double> shape(nloc);
    double io = inside_element(id, nodes, shape);
    
    if(nearest.first==-1){
      nearest.first = id;
      nearest.second = io;
      for(size_t i=0;i<nloc;i++)
        nearest_shape[i] = shape[i];
    }else{
      if(io>nearest.second){
        nearest.first = id;
        nearest.second = io;
        for(size_t i=0;i<nloc;i++)
          nearest_shape[i] = shape[i];
      }
    }
    
    delete [] nodes;
  }
  
  void visitData(std::vector<const SpatialIndex::IData*>& v) {}
  
private:
  double inside_element(int id, const int *n, vector<double>& shape){
    if((element_type==TRI)&&(dimension==2))
      return inside_triangle_2d(id, n, shape);
    else if(element_type==TETRA)
      return inside_tetra(id, n, shape);
    // else
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): Gerard was regrettably kidnapped "
        <<"by aliens before he found time to write support for element type "<<element_type<<endl;
    exit(-1);      
    
    // shup up compiler warning
    return 0;
  }
  
  double inside_triangle_2d(int id, const int *n, vector<double>& shape){
    double Normal = normal2d(sx[n[0]], sy[n[0]], sx[n[1]], sy[n[1]], sx[n[2]], sy[n[2]]);
    
    double N=normal2d(sx[n[0]], sy[n[0]],
                      sx[n[1]], sy[n[1]],
                      find_x, find_y);
    shape[2] = N/Normal;
    
    N = normal2d(sx[n[1]], sy[n[1]],
                 sx[n[2]], sy[n[2]],
                 find_x, find_y);
    shape[0] = N/Normal;
    
    N = normal2d(sx[n[2]], sy[n[2]],
                 sx[n[0]], sy[n[0]],
                 find_x, find_y);
    shape[1] = N/Normal;
    
    return min(shape[0], min(shape[1], shape[2]));
  }

  double inside_tetra(int id, const int *n, vector<double>& shape){
    double det = calc_vol(sx[n[0]], sy[n[0]], sz[n[0]],
                          sx[n[1]], sy[n[1]], sz[n[1]],
                          sx[n[2]], sy[n[2]], sz[n[2]],
                          sx[n[3]], sy[n[3]], sz[n[3]]);
    
    shape[0] = calc_vol(sx[n[1]], sy[n[1]], sz[n[1]],
                        sx[n[3]], sy[n[3]], sz[n[3]],
                        sx[n[2]], sy[n[2]], sz[n[2]],
                        find_x, find_y, find_z)/det;
    
    shape[1] = calc_vol(sx[n[0]], sy[n[0]], sz[n[0]],
                        sx[n[2]], sy[n[2]], sz[n[2]],
                        sx[n[3]], sy[n[3]], sz[n[3]],
                        find_x, find_y, find_z)/det;
    
    shape[2] = calc_vol(sx[n[0]], sy[n[0]], sz[n[0]],
                        sx[n[3]], sy[n[3]], sz[n[3]],
                        sx[n[1]], sy[n[1]], sz[n[1]],
                        find_x, find_y, find_z)/det;
    
    shape[3] = calc_vol(sx[n[0]], sy[n[0]], sz[n[0]],
                        sx[n[1]], sy[n[1]], sz[n[1]],
                        sx[n[2]], sy[n[2]], sz[n[2]],
                        find_x, find_y, find_z)/det;
    
    return min(min(shape[0], shape[1]), min(shape[2], shape[3]));
  }
  
  double normal2d(double x0, double y0,
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

  double calc_vol(double x0, double y0, double z0,
                  double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double x3, double y3, double z3) const{
    
    double a1 = x1-x0;
    double a2 = y1-y0;
    double a3 = z1-z0;
    
    double b1 = x2-x0;
    double b2 = y2-y0;
    double b3 = z2-z0;
    
    double c1 = x3-x0;
    double c2 = y3-y0;
    double c3 = z3-z0;
    
    // volume = | r_a r_b r_c | / 6
    return a1*(b2*c3 - b3*c2) - b1*(a2*c3 - a3*c2) + c1*(a2*b3 - a3*b2);
  }
  
  const double *sx, *sy, *sz;
  unsigned long element_type, dimension, nloc;
  pair<int, double> nearest;
  vector<double> nearest_shape;
  double find_x, find_y, find_z;

  const unsigned long TRI, QUAD, TETRA, HEX;
};

Picker::Picker(const double& tol): TRI(5), QUAD(9), TETRA(10), HEX(12){
#ifdef HAVE_VTK
  assert(TRI==VTK_TRIANGLE);
  assert(QUAD==VTK_QUAD);
  assert(TETRA==VTK_TETRA);
  assert(HEX==VTK_HEXAHEDRON);
#endif

  verbose = false;
  spherical_geometry = false;
  dimension = 0;
  edimension = 0;
  
  rtree_north = NULL;
  rtree_south = NULL;
  rtree = NULL;

  storage_manager_north = NULL;
  storage_manager_south = NULL;
  storage_manager = NULL;

  indexId_north=1;
  indexId_south=2;
  indexId=3;

  // R-Tree parameters - need to see how these need to be optimized
  // for 2D and 3D.
  fillFactor = 0.7;
  indexCapacity = 100;
  leafCapacity = 100;
  variant = SpatialIndex::RTree::RV_RSTAR;
  
  this->SetOwnershipTolerance(tol);
}

Picker::~Picker(){
  if(verbose)
    cout<<"Picker::~Picker()\n";

  if(rtree_north)
    delete rtree_north;
  if(rtree_south)
    delete rtree_south;
  if(rtree)
    delete rtree;

  if(storage_manager_north)
    delete storage_manager_north;
  if(storage_manager_south)
    delete storage_manager_south;
  if(storage_manager)
    delete storage_manager;
}

void Picker::CreateStorageManagers(){
  if(verbose)
    cout<<"void Picker::CreateStorageManagers()\n";

  if(dimension==0){
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): "
        <<"Source data set has not been established\n";
    exit(-1);
  }
  
  data_block_size = nloc*sizeof(int);

  if(spherical_geometry){
    if(element_type==TETRA)
      edimension = 3;
    else
      edimension = 2;
  }else{
    edimension = dimension;
  }

  // Create storage managers
  if(spherical_geometry){
    storage_manager_north = SpatialIndex::StorageManager::createNewMemoryStorageManager();
    storage_manager_south = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  }else{
    storage_manager = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  }
}

void Picker::DumpVTU(const char *basename) const{
  if(verbose)
    cout<<"void Picker::DumpVTU("<<basename<<") const\n";
  
  if(dimension==0)
    return;

#ifdef HAVE_VTK
  if(elements_north.size()>0){ // north
    string filename=string(basename)+"_north.vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();
    
    size_t ecnt = element_ids_north.size();
    size_t ncnt = sx.size();

    // Point definitions
    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    for(size_t i=0; i<ncnt; i++){
      if(edimension==2)
        newPts->InsertNextPoint(sx[i], sy[i], 0.0);
      else
        newPts->InsertNextPoint(sx[i], sy[i], sz[i]);
    }
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();
    
    // Cell definitions
    vtkIdType Cell[nloc];
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<nloc; j++){
        assert(elements_north[i*nloc+j]>=0);
        assert(elements_north[i*nloc+j]<(int)ncnt);
        Cell[j] = elements_north[i*nloc+j];
      }
      ug->InsertNextCell(element_type, nloc, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }

  if(elements_south.size()>0){ // south
    string filename=string(basename)+"_south.vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();

    const size_t ecnt = element_ids_north.size();

    // Point definitions
    size_t ncnt = 0;
    for(size_t i=0; i<ecnt; i++)
      for(size_t j=0; j<nloc; j++)
        ncnt = max(ncnt, (size_t)elements_south[i*nloc+j]);
    ncnt++;

    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    for(size_t i=0; i<ncnt; i++)
      if(edimension==2)
        newPts->InsertNextPoint(sx[i], sy[i], 0.0);
      else
        newPts->InsertNextPoint(sx[i], sy[i], sz[i]);
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();

    // Cell definitions
    vtkIdType Cell[nloc];
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<nloc; j++){
        assert(elements_south[i*nloc+j]>=0);
        assert(elements_south[i*nloc+j]<(int)ecnt);
        Cell[j] = elements_south[i*nloc+j];
      }
      ug->InsertNextCell(element_type, nloc, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }

  if(elements.size()>0){ // plane
    string filename=string(basename)+".vtu";
    if(verbose)
      cout<<"Dumping out "<<filename<<endl;

    vtkUnstructuredGrid *ug=vtkUnstructuredGrid::New();
    
    const size_t ecnt = elements.size()/nloc;
    
    // Point definitions
    size_t ncnt = 0;
    for(size_t i=0; i<ecnt; i++)
      for(size_t j=0; j<nloc; j++)
        ncnt = max(ncnt, (size_t)elements[i*nloc+j]);
    ncnt++;

    vtkPoints *newPts = vtkPoints::New();
    newPts->SetDataTypeToDouble();
    for(size_t i=0; i<ncnt; i++)
      if(dimension==2)
        newPts->InsertNextPoint(_x[i], _y[i], 0.0);
      else
        newPts->InsertNextPoint(_x[i], _y[i], _z[i]);
    ug->SetPoints(newPts);
    ug->Modified();
    ug->Update();
    newPts->Delete();
    
    // Cell definitions
    vtkIdType Cell[ncnt];
    for(size_t i=0; i<ecnt; i++){
      for(size_t j=0; j<nloc; j++){
        assert(elements[i*nloc+j]>=0);
        assert(elements[i*nloc+j]<(int)ecnt);
        Cell[j] = elements[i*nloc+j];
      }
      ug->InsertNextCell(element_type, nloc, Cell);
    }
    
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(ug);
    writer->Write();
  }
#endif
}

void Picker::SetOwnershipTolerance(const double& tol){
  assert(tol >= 0.0);
  this->tol = tol;

  return;
}

double Picker::GetOwnershipTolerance() const{
  return tol;
}

// Returns 0 if an element containing the queried point is found,
// otherwise -1, indicating that the nearest element is being passed
// back.
int Picker::Find(double x, double y, int *eid, double *shape){
  if(verbose)
    cout<<"int Picker::Find("<<x<<", "<<y<<", int *eid, double *shape)\n";

  assert(edimension==2);
  double r[2];
  r[0] = x; r[1] = y;
  return Find(r, eid, shape);
}

int Picker::Find(double x, double y, double z, int *eid, double *shape){
  if(verbose)
    cout<<"int Picker::Find("<<x<<", "<<y<<", "<<z<<", int *eid, double *shape)\n";

  assert(edimension==3);

  double r[3];
  r[0] = x; r[1] = y; r[2] = z;
  return Find(r, eid, shape);
}

int Picker::Find(const double *x, int *eid, double *shape){
  if(verbose)
    cout<<"int Picker::Find(const double *x, int *eid, double *shape)\n";

  *eid = -1;
  double pos[3];
  
  if(spherical_geometry){
    assert(dimension==3);
    if(x[2]>=0){
      projection::cartesian2stereographic(0.0, 90.0, x[0], x[1], x[2], pos[0], pos[1], pos[2]);
      if(verbose)
        cout<<"Searching northern hemisphere for "<<pos[0]<<", "<<pos[1]<<", "<<pos[2]<<endl;
    }else{
      projection::cartesian2stereographic(0.0,-90.0, x[0], x[1], x[2], pos[0], pos[1], pos[2]);
      if(verbose)
        cout<<"Searching southern hemisphere for "<<pos[0]<<", "<<pos[1]<<", "<<pos[2]<<endl;
    }
  }else{
    pos[0] = x[0];
    pos[1] = x[1];
    pos[2] = x[2];
  }
  
  SpatialIndex::Point p = SpatialIndex::Point(pos, edimension);  
  SIVisitor v(element_type, edimension);
    
  if(spherical_geometry){
    if(edimension==2)
      v.SetCoordinates(&(sx[0]), &(sy[0]));
    else
      v.SetCoordinates(&(sx[0]), &(sy[0]), &(sz[0]));
  }else{
    if(dimension==2)
      v.SetCoordinates(&(_x[0]), &(_y[0]));
    else
      v.SetCoordinates(&(_x[0]), &(_y[0]), &(_z[0]));
  }

  v.SetPoint(pos);
    
  if(spherical_geometry){
    if(x[2]>=0){
      rtree_north->pointLocationQuery(p, v);
    }else{
      rtree_south->pointLocationQuery(p, v);
    }
  }else{
    rtree->pointLocationQuery(p, v);
  }

  *eid = v.GetHit(shape);
  
  double min_shape=shape[0];
  for(unsigned long i=1; i<nloc; i++)
    min_shape = min(min_shape, shape[i]);
  
  if(verbose){
    cout<<"Hit element = "<<*eid<<"\n"
          "Minimum interpolation weight = "<<min_shape<<endl;
  }
  
  if(min_shape<-tol)
    return -1;
  
  return 0;
}

bool Picker::IsParallel() const{
  if(verbose)
    cout<<"bool Picker::IsParallel() const\n";
  
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
int Picker::pFind(const double *x, const double *y, const double *z, int nnodes,
                  int *eid, double *shape, int *host, int different_domains){
  if(verbose)
    cout<<"int Picker::pFind(const double *x, const double *y, const double *z, int nnodes, int *eid, double *shape, int *host)\n";

  size_t myrank = 0;
  size_t NProcs = 1;
  
#ifdef HAVE_MPI
  if(IsParallel()){
    myrank = MPI::COMM_WORLD.Get_rank();
    NProcs = MPI::COMM_WORLD.Get_size();
  }
#endif
  
  vector<double> missing_xyz;
  vector<int> lut;
  for(size_t i=0;i<(size_t)nnodes;i++){
    int stat;
    if(z){
      if(verbose)
        cout<<"searching for node "<<i<<"(3d): "<<x[i]<<", "<<y[i]<<", "<<z[i]<<endl;
      stat = Find(x[i], y[i], z[i], eid+i, shape+i*nloc);
    }else{
      if(verbose)
        cout<<"searching for node "<<i<<"(2d): "<<x[i]<<", "<<y[i]<<", "<<endl;
      stat = Find(x[i], y[i], eid+i, shape+i*nloc);
    }

    if(verbose){
      if(stat){
        cout<<"node not found\n";
      }else{
        cout<<"resulting shape: ";
        for(size_t j=0;j<nloc;j++)
          cout<<shape[i*nloc+j]<<" ";
        cout<<endl;
      }
    }

    double min_weight = shape[i*nloc];
    for(size_t j=1;j<nloc;j++)
      min_weight = min(min_weight, shape[i*nloc+j]);

    if(min_weight<-tol){
      missing_xyz.push_back(x[i]);
      missing_xyz.push_back(y[i]);
      if(z)
        missing_xyz.push_back(z[i]);
      lut.push_back(i);
      host[i] = -1;
    }else{
      host[i] = myrank;
    }
  }

#ifdef HAVE_MPI
  if(IsParallel() && (different_domains != 1)){
    for(size_t i=0;i<NProcs;i++){
      int num_missing=lut.size();
      MPI::COMM_WORLD.Bcast(&num_missing, 1, MPI::INT, i);
      
      if(num_missing==0)
        continue;

      vector<double> looking_for_xyz;
      looking_for_xyz.resize(num_missing*edimension);
      if(myrank==i)
        looking_for_xyz = missing_xyz;
      
      MPI::COMM_WORLD.Bcast(&(looking_for_xyz[0]), num_missing*edimension, MPI::DOUBLE, i);
    
      if(myrank!=i){
        size_t missing_count = looking_for_xyz.size()/edimension;
        vector<int> found; // query point number
        for(size_t j=0;j<missing_count;j++){
          vector<double> missing_shape(nloc), found_shape;
          int missing_eid;
          Find(&(looking_for_xyz[j*edimension]), &missing_eid, &(missing_shape[0]));

          double min_shape=shape[0];
          for(unsigned long i=1; i<nloc; i++)
            min_shape = min(min_shape, shape[i]);
          
          if(min_shape<0.0)
            found.push_back(j);
        }
        MPI::COMM_WORLD.Send(&(found[0]), found.size(), MPI::INT, i, 0);
      }else{
        for(size_t p=0;p<NProcs-1;p++){
          MPI::Status status;
          MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, 0, status);
          int rcnt = status.Get_count(MPI::INT);
          int source = status.Get_source();
          vector<int> buff(rcnt);
          MPI::COMM_WORLD.Recv(&(buff[0]), rcnt, MPI::INT, source, 0);
          for(vector<int>::const_iterator it=buff.begin();it!=buff.end();++it){
            // Perhaps it's already been found by another host
            if(host[lut[*it]]<0)
              host[lut[*it]] = source;
          }
        }
      }
    }
  }
#endif
  
  int err=0;
  for(size_t i=0;i<(size_t)nnodes;i++)
    if(host[i]<0){
      if(verbose){
        if(z){
          cerr<<"For node "<<i<<"(3d): "<<x[i]<<", "<<y[i]<<", "<<z[i]<<endl;
        }else{
          cerr<<"For node "<<i<<"(2d): "<<x[i]<<", "<<y[i]<<", "<<endl;
        }
        cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__
            <<"): point not contained within mesh.\n";
      }
      err = -1;
    }
  
  return err;
}

void Picker::SetSource(const double *x, const double *y,
                       const int *enlist, int nelements,
                       int _element_type, int zero_index){
  if(verbose)
    cout<<"void Picker::SetSource(const double *x, const double *y, const int *enlist, int nelements, int _element_type, int zero_index)\n";
  
  SetSource(x, y, NULL,
            enlist, nelements,
            _element_type, zero_index);
}

void Picker::SetSource(const double *x, const double *y, const double *z,
                       const int *enlist, int nelements,
                       int _element_type, int zero_index){
  if(verbose)
    cout<<"void Picker::SetSource(const double *x, const double *y, const double *z, const int *enlist, int nelements, int _element_type, int zero_index)\n";
  
  if(z==NULL)
    dimension = 2;
  else
    dimension = 3;
  
  element_type = _element_type;
  if(verbose)
    cout<<"Element type: ";
  if(element_type==TRI){
    if(verbose)
      cout<<"TRI\n";
    nloc = 3;
  }else if(element_type==QUAD){
    if(verbose)
      cout<<"QUAD\n";
    nloc = 4;
  }else if(element_type==TETRA){
    if(verbose)
      cout<<"TETRA\n";
    nloc = 4;
  }else if(element_type==HEX){
    if(verbose)
      cout<<"HEX\n";
    nloc = 8;
  }else{
    cerr<<"ERROR ("<<__FILE__<<", "<<__LINE__<<"): unknown element type in Picker.\n";
    exit(-1);
  }

  CreateStorageManagers();
 
  if(spherical_geometry){
    // Create two sets of elements - those that has nodes in the
    // northern helisphere and those with elements in the southern
    // hemisphere. Note that elements may appear in both sets and their
    // nodes will be duplicated. Each hemisphere will be mapped with a
    // stereographic projection using the opposing pole as the
    // projection point which goes to infinity.
    
    set<int> north_elm_set, south_elm_set;
    map<int, int> renumber_north, renumber_south;
    int mxnid=0;
    for(size_t i=0;i<(size_t)nelements;i++){
      bool north=false, south=false;
      
      for(size_t j=0;j<nloc;j++){
        int nid = enlist[i*nloc+j]-zero_index;
        if(z[nid]>=0)
          north = true;
        else
          south = true;
      }
      
      if(north){
        north_elm_set.insert(i);
        for(size_t j=0;j<nloc;j++){
          int nid = enlist[i*nloc+j]-zero_index;
          renumber_north[nid] = -1;
          mxnid=max(nid, mxnid);
        }
      }
      
      if(south){    
        south_elm_set.insert(i);
        for(size_t j=0;j<nloc;j++){
          int nid = enlist[i*nloc+j]-zero_index;
          renumber_south[nid] = -1;
          mxnid=max(nid, mxnid);
        }
      }
    }
        
    // Finalise renumbering
    int pos=0;
    for(map<int, int>::iterator it=renumber_north.begin(); it!=renumber_north.end(); ++it){
      it->second = pos++;
    }
    for(map<int, int>::iterator it=renumber_south.begin(); it!=renumber_south.end(); ++it){
      it->second = pos++;
    }
    
    // Store renumbered elements
    for(set<int>::const_iterator it=north_elm_set.begin(); it!=north_elm_set.end(); ++it){
      element_ids_north.push_back(*it);
      for(size_t j=0;j<nloc;j++){
        int nid = enlist[(*it)*nloc+j]-zero_index;
        elements_north.push_back(renumber_north[nid]);
      }
    }
    
    for(set<int>::const_iterator it=south_elm_set.begin(); it!=south_elm_set.end(); ++it){
      element_ids_south.push_back(*it);
      for(size_t j=0;j<nloc;j++){
        int nid = enlist[(*it)*nloc+j]-zero_index;
        elements_south.push_back(renumber_south[nid]);
      }
    }
    
    // Project and store all the valid points
    if(element_type==TRI){
      double xx, yy;
      for(map<int, int>::const_iterator it=renumber_north.begin(); it!=renumber_north.end(); ++it){
        projection::cartesian2stereographic(0.0, 90.0, x[it->first], y[it->first], z[it->first], xx, yy);
        sx.push_back(xx);
        sy.push_back(yy);
      }
      
      for(map<int, int>::const_iterator it=renumber_south.begin(); it!=renumber_south.end(); ++it){
        projection::cartesian2stereographic(0.0, -90.0, x[it->first], y[it->first], z[it->first], xx, yy);
        sx.push_back(xx);
        sy.push_back(yy);
      }
    }else{
      double xx, yy, zz;
      for(map<int, int>::const_iterator it=renumber_north.begin(); it!=renumber_north.end(); ++it){
        projection::cartesian2stereographic(0.0, 90.0, x[it->first], y[it->first], z[it->first], xx, yy, zz);
        sx.push_back(xx);
        sy.push_back(yy);
        sz.push_back(zz);
      }
      
      for(map<int, int>::const_iterator it=renumber_south.begin(); it!=renumber_south.end(); ++it){
        projection::cartesian2stereographic(0.0, -90.0, x[it->first], y[it->first], z[it->first], xx, yy, zz);
        sx.push_back(xx);
        sy.push_back(yy);
        sz.push_back(zz);
      }
    }
    
    // Fill RTree's
    if(element_type==TRI){
      PickerMeshDataStream nstream(&(sx[0]), &(sy[0]), &(elements_north[0]), 0,
                                   element_ids_north.size(), nloc, &(element_ids_north[0]));
      rtree_north = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, nstream, *storage_manager_north, fillFactor,
                                                                   indexCapacity, leafCapacity, 2, variant, indexId_north);        
      
      PickerMeshDataStream sstream(&(sx[0]), &(sy[0]), &(elements_south[0]), 0,
                                   element_ids_south.size(), nloc, &(element_ids_south[0]));
      rtree_south = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, sstream, *storage_manager_south, fillFactor,
                                                                   indexCapacity, leafCapacity, 2, variant, indexId_south);
    }else{
      PickerMeshDataStream nstream(&(sx[0]), &(sy[0]), &(sz[0]), &(elements_north[0]), 0,
                                   element_ids_north.size(), nloc, &(element_ids_north[0]));
      rtree_north = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, nstream, *storage_manager_north, fillFactor,
                                                                   indexCapacity, leafCapacity, 3, variant, indexId_north);
      
      PickerMeshDataStream sstream(&(sx[0]), &(sy[0]), &(sz[0]), &(elements_south[0]), 0, 
                                   element_ids_south.size(), nloc, &(element_ids_south[0]));
      rtree_south = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, sstream, *storage_manager_south, fillFactor,
                                                                   indexCapacity, leafCapacity, 3, variant, indexId_south);
    }
  }else{
    _x = x;
    _y = y;
    _z = z;
    
    // Fill RTree
    PickerMeshDataStream stream(x, y, z, enlist, 1, nelements, nloc);
    rtree = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, stream, *storage_manager, fillFactor,
                                                           indexCapacity, leafCapacity, edimension, variant, indexId);
  }
}

void Picker::SetSphericalGeometry(bool stat){
  if(verbose)
    cout<<"void Picker::SetSphericalGeometry(bool stat)\n";

  spherical_geometry = stat;
}

map<int, Picker> pickers;

extern "C"{
#define picker_create_fc F77_FUNC(picker_create, PICKER_CREATE)
#define picker_destroy_fc F77_FUNC(picker_destroy, PICKER_destroy)
#define picker_find2d_fc F77_FUNC(picker_find2d, PICKER_FIND2D)
#define picker_find3d_fc F77_FUNC(picker_find3d, PICKER_FIND3D)
#define picker_pfind2d_fc F77_FUNC(picker_pfind2d, PICKER_PFIND2D)
#define picker_pfind3d_fc F77_FUNC(picker_pfind3d, PICKER_PFIND3D)
  
  void picker_create_fc(const double x[], const double y[], const double z[],
                        const int enlist[], const int *nelm, const int *element_type,
                        const int *dimension, const int *fortran_zero_index, const int *flat_earth,
                        int *picker_id){
    *picker_id = 1;
    while(pickers.find(*picker_id)!=pickers.end())
      *picker_id = *picker_id + 1;
    
    pickers[*picker_id].SetSphericalGeometry(*flat_earth==0);
    
    if(*dimension==2)
      pickers[*picker_id].SetSource(x, y, NULL, 
                                    enlist, *nelm,
                                    *element_type, *fortran_zero_index);
    else
      pickers[*picker_id].SetSource(x, y, z, 
                                    enlist, *nelm,
                                    *element_type, *fortran_zero_index);
    return;
  }

  void picker_destroy_fc(const int *picker_id){
    pickers.erase(*picker_id);
  }

  void picker_find2d_fc(const int *picker_id, const double *x, const double *y,
                        int *eid, double shape[], const double *tol){
    Picker *picker = &pickers[*picker_id];
    
    picker->SetOwnershipTolerance(*tol);                   
    picker->Find(*x, *y, eid, shape);
    *eid = *eid + 1;
    return;
  }

  void picker_find3d_fc(const int *picker_id, const double *x, const double *y, const double *z,
                        int *eid, double shape[], const double *tol){
    Picker *picker = &pickers[*picker_id];
    
    picker->SetOwnershipTolerance(*tol);  
    picker->Find(*x, *y, *z, eid, shape);
    *eid = *eid + 1;
    return;
  }

  void picker_pfind2d_fc(const int *picker_id, const double x[], const double y[], const int *nnodes,
                         int eid[], double shape[], int rank[], const double *tol, const int *different_domains){
    Picker *picker = &pickers[*picker_id];
    
    picker->SetOwnershipTolerance(*tol);              
    picker->pFind(x, y, NULL, *nnodes, eid, shape, rank, *different_domains);
    for(int i=0;i<*nnodes;i++)
      eid[i]++;
    return;
  }

  void picker_pfind3d_fc(const int *picker_id, const double x[], const double y[], const double z[], const int *nnodes,
                         int eid[], double shape[], int rank[], const double *tol, const int *different_domains){
    Picker *picker = &pickers[*picker_id];
    
    picker->SetOwnershipTolerance(*tol); 
    picker->pFind(x, y, z, *nnodes, eid, shape, rank, *different_domains);
    for(int i=0;i<*nnodes;i++)
      eid[i]++;    
    return;
  }

}

void Picker::VerboseOff(){
  verbose = false;
}

void Picker::VerboseOn(){
  verbose = true;
}

#ifdef TEST_PICKER

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
  
  Picker vertical;
  vertical.VerboseOn();
  vertical.SetSphericalGeometry(true);
  vertical.SetSource(&(x[0]), &(y[0]), &(z[0]), &(ENList[0]), ecnt, 5, 0);  
  vertical.DumpVTU("projection");

  vector<double> shape(ncnt*3);
  vector<int> eid(ncnt), host(ncnt);
  vertical.pFind(&(x[0]), &(y[0]), &(z[0]), 3,
                 &(eid[0]), &(shape[0]), &(host[0]), 0);

  for(vtkIdType i=0;i<3;i++){
    cout<<"Looking for element containing node: "<<i
        <<", ("<<x[i]<<", "<<y[i]<<", "<<z[i]<<")"<<endl;
    
    cout<<"Found element: "<<eid[i]<<". Contains:\n";
    for(size_t j=0;j<3;j++){
      int nid = ug->GetCell(eid[i])->GetPointId(j);
      cout<<nid<<":\t"<<x[nid]<<" "<<y[nid]<<" "<<z[nid]<<endl;
    }
  }
  
#else
  cerr<<"ERROR: no vtk support compiled\n";
#endif
  
  return 0;
}
#endif
