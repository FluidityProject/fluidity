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

#include "Adaptivity.h"
#include "cinterfaces.h"

using namespace std;

bool Adaptivity::verbose=false;

// Constructer
Adaptivity::Adaptivity(){
  if(verbose)
    cout<<"Adaptivity::Adaptivity()\n";

  newNNodes     = 0;
  newNElements  = 0;
  newNSElements = 0;

  newX = NULL;
  newY = NULL;
  newZ = NULL;

  numberingOffset = 1;
  
  volumeID = NULL;

  NPrivateNodes = 0;
  NProcs  = getNProcessors();
  Gather  = NULL;
  Scatter = NULL;
  NGather = 0;
  NScatter= 0;
  ATOSEN  = NULL;
  ATOREC  = NULL;

  set_adapt_sweeps(10);

  AdaptOpts[1] = LOGICAL_TRUE; // collapse edges if true
  disableGeometryDiscovery();
  enableMovementConnectivity();
  enableRefinement();
  disableSurfaceLock();

  setFunctionalTolerance(0.7);

  nloc  = 4;
  snloc = 3;
  
  NWSZEN = 0;
  NWSZSN = 0;
  NWSZNN = 0;
  NWNDLC = 0;
  NWSROW = 0;
  NWENLB = 0;
  NWENLS = 0;
  NWSNLB = 0; 
  NWSNLS = 0; 
  NWSFID = 0;
  NWELRG = 0;
  NWNODX = 0; 
  NWNODY = 0;
  NWNODZ = 0;
  NEWMTX = 0;
  NEWFLD = 0;
  ADPBIG = 0; 
  ADPNOD = 0;

  interpolate = false;
};

afloat_t Adaptivity::edgeLengthDistribution(afloat_t w){
  if(verbose)
    cout<<"float Adaptivity::edgeLengthDistribution(" << w << ")" << endl;
#ifdef HAVE_VTK
  int edgecount = 0;
  int badedges = 0;
  const int dim = 3;
  
  // Figure out the edges.
  std::set<std::pair<int, int> > edges;

  // Form NNlist.
  for(int e=0;e<ug->GetNumberOfCells();e++){
    int k, l;
    vtkTetra *tetra = (vtkTetra *)ug->GetCell(e);
    for(size_t i=0;i<4;i++){
      for(size_t j=i+1;j<4;j++){
        k = tetra->GetPointId(i);
        l = tetra->GetPointId(j);
        edges.insert(std::pair<int, int>(min(k, l), max(k, l)));
      }
    }
  }
  
  for (std::set<std::pair<int, int> >::iterator edge = edges.begin(); edge != edges.end(); edge++){
    edgecount++;
    // First, compute the metric over the edge by averaging the nodal metrics.
    afloat_t edge_metric[dim*dim];
    for (int i = 0; i < dim; i++){
      for (int j = 0; j < dim; j++){
        edge_metric[i*dim+j] = (Metric[edge->first*dim*dim + i*dim + j] + Metric[edge->second*dim*dim + i*dim + j]) / 2.0;
      }
    }
    
    // Now, compute the inner product.
    afloat_t vec[dim];
    vec[0] = X[edge->first] - X[edge->second];
    vec[1] = Y[edge->first] - Y[edge->second];
    vec[2] = Z[edge->first] - Z[edge->second];

    afloat_t len = 0;
    for (int i = 0; i < dim; i++){
      afloat_t tmp = 0;
      for (int j = 0; j < dim; j++)
        tmp = tmp + edge_metric[i*dim + j] * vec[j];
      len = len + vec[i] * tmp;
    }
    len = sqrt(len);
    if ((len < 1.0 - w) || (len > 1.0 + w))
      badedges++;
  }
  
  return badedges / (afloat_t)edgecount;
#else
  cerr<<"WARNING: VTK support has not been compiled.\n";
  return 0.0;
#endif
}

void Adaptivity::adapt(){
  if(verbose)
    cout<<"void Adaptivity::adapt()\n";

  //
  // Hardwired stuff
  //
  
  // Only does 3D
  int Geom3D = LOGICAL_TRUE;
  
  // DOTOP is the minimum element functional to allow to be changed (min 0.1)
  afloat_t DOTOP  = (fabs(MESTP1)>0.15)?fabs(MESTP1):0.15;

  //  MINCHG is the minimum relative change allowed to a local functional (min 0.01)
  afloat_t MINCHG = 0.01;

  // if true then use the quality of the element rather than insphere
  int USEQ = LOGICAL_FALSE;
    
  NWSZEN = 0;
  NWSZSN = 0;
  NWSZEN = 0;
  NWSZSN = 0;
  
  // These variables...if they are set to a non-negitive number they
  // write back the data - otherwise no
  NEWMTX =  1; // New metric values
  NWNDLC =  1; // Indicates which old element contains which new node
  NWSROW = -1;
  
  //
  // Local variables
  //
  
  // Find the expected number of elements given a metric
  int NNodes = X.size();
  int NElements = ENList.size()/nloc;
  int NSElements = SENList.size()/snloc;
  int XPCTEL = get_predicted_nelements_fc(&(Metric[0]), &(X[0]), &(Y[0]), &(Z[0]),
                                          &(ENList[0]), &NNodes, &NElements, &nloc);

  // Calculate how much memory adaptivity needs and allocate it.
  int intBuffer_size, floatBuffer_size;
  int DUM1=-1, DUM2=-1,DUM3=-1,one=1;
  int flag=LOGICAL_TRUE;
  adaptmem_(&NNodes, &NElements, &DUM1, &NSElements, &one,
	    &XPCTEL, &DUM2, &DUM3, &flag, &intBuffer_size, &floatBuffer_size);

  intBuffer.resize(intBuffer_size);
  floatBuffer.resize(floatBuffer_size);
  
  vector<int> ENLBasePtr(NElements+1);
  for(int i=0; i<=NElements; i++)
    ENLBasePtr[i] = i*nloc;
  int sizeENList = ENList.size();
  
  vector<int> SENLBasePtr(NSElements+1);
  for(int i=0; i<=NSElements; i++)
    SENLBasePtr[i]=i*snloc;
  int sizeSENList = SENList.size();
  
  vector<int> ElementRegion;
  const int *_volumeID = volumeID;
  if(volumeID==NULL){
    ElementRegion.resize(NElements);
    _volumeID=&(ElementRegion[0]);
    for(int i=0; i<NElements; i++)
      ElementRegion[i] = 1;
  }
  
  if(surfID.size()==0){
    // Probably should calculate this using the co-planar patch
    // algorithm if it's not already done.
    surfID.resize(NSElements);
    for(int i=0;i<NSElements;i++)
      surfID[i]=1; 
  }
  
  // Halo information gets over written so...
  int *_newGather=NULL, *_newScatter=NULL, *_newATOSEN=NULL, *_newATOREC=NULL;
  if(NProcs>1){
    newGather.resize(NGather);
    _newGather = &(newGather[0]);
    memcpy(_newGather, Gather, NGather*sizeof(int));
    
    newScatter.resize(NScatter);
    _newScatter = &(newScatter[0]);
    memcpy(_newScatter, Scatter, NScatter*sizeof(int));
    
    newATOSEN.resize(NProcs+1);
    _newATOSEN = &(newATOSEN[0]);
    memcpy(_newATOSEN, ATOSEN, (NProcs+1)*sizeof(int));
    
    newATOREC.resize(NProcs+1);
    _newATOREC = &(newATOREC[0]);
    memcpy(_newATOREC, ATOREC, (NProcs+1)*sizeof(int));
  }

  newNPrivateNodes = NPrivateNodes;

  int debug_level = 0;
  int *PRDNDS = NULL, NPRDND=0;

  //if(verbose){
  //  chcnsy = LOGICAL_TRUE;
  //  dbg    = LOGICAL_TRUE;
  //}else{
  chcnsy = LOGICAL_FALSE;
  dbg    = LOGICAL_FALSE;
  //}
  
  // Gerard: this appears to be bugged
  // -1 lets libadapt figure it out for itself
  //int mxnods = 2.0*XPCTEL * (NNodes / (double)NElements);
  int mxnods = -1;

  int totfre=0;
  for(vector<int>::const_iterator it=nfreedom.begin();it!=nfreedom.end();++it){
    totfre+=*it;
  }

  // Call the main adaptivity routine  
  adptvy_(&(intBuffer[0]), &intBuffer_size, &(floatBuffer[0]), &floatBuffer_size,
	  &Geom3D, &SRFGMY, &USEQ,
	  &NNodes, &NElements, &NSElements, &mxnods,
	  &sizeENList, &(ENLBasePtr[0]), &(ENList[0]), _volumeID,
	  &CLCGMY, 
	  &sizeSENList, &(SENLBasePtr[0]), &(SENList[0]), &(surfID[0]),
	  PRDNDS, &NPRDND,
	  &(X[0]), &(Y[0]), &(Z[0]),
	  &NNodes, &NElements, &sizeENList, &(ENList[0]), &(ENLBasePtr[0]),
	  &(X[0]), &(Y[0]), &(Z[0]),
	  &(Metric[0]), &(fields[0]), &(nfreedom[0]), &totfre, &nfields,
	  &XPCTEL,
	  &newNNodes, &newNElements, &newNSElements,
	  &NWSZEN, &NWSZSN, &NWSZNN, &NWNDLC, &NWSROW,
	  &NWENLB, &NWENLS, &NWSNLB, &NWSNLS, &NWSFID,
	  &NWELRG, &NWNODX, &NWNODY, &NWNODZ,
	  &NEWMTX, &NEWFLD,
	  &ADPBIG, &ADPNOD,
	  &DOTOP,  &MINCHG, &MaxNumberAdaptIterations, AdaptOpts, &TWOSTG, &TOGTHR,
	  _newGather, _newScatter, &NGather, &NScatter, &newNPrivateNodes,
	  _newATOSEN, _newATOREC,  &NProcs, &debug_level, &dbg, &chcnsy);

  // Zero tolerence for errors from adaptivity
  assert(newNNodes>0);
  
  newX = &(floatBuffer[NWNODX-1]);
  newY = &(floatBuffer[NWNODY-1]);
  newZ = &(floatBuffer[NWNODZ-1]);
}

void Adaptivity::enableGeometryDiscovery(){
  if(verbose)
    cout<<"void Adaptivity::enableGeometryDiscovery()\n";
  CLCGMY = LOGICAL_TRUE;
}

void Adaptivity::enableMovementConnectivity(){
  if(verbose)
    cout<<"void Adaptivity::enableMovementConnectivity()\n";
  TOGTHR = LOGICAL_TRUE;
}

void Adaptivity::enableRefinement(){
  if(verbose)
    cout<<"void Adaptivity::enableRefinement()\n";
  TWOSTG = LOGICAL_FALSE;

  AdaptOpts[0] = LOGICAL_TRUE;
  AdaptOpts[2] = LOGICAL_TRUE;
  AdaptOpts[3] = LOGICAL_TRUE;
  AdaptOpts[4] = LOGICAL_TRUE;
  AdaptOpts[5] = LOGICAL_TRUE;
}

void Adaptivity::enableSurfaceLock(){
  if(verbose)
    cout<<"void Adaptivity::enableSurfaceLock()\n";
  SRFGMY = LOGICAL_TRUE;
}

void Adaptivity::disableGeometryDiscovery(){
  if(verbose)
    cout<<"void Adaptivity::disableGeometryDiscovery()\n";
  CLCGMY = LOGICAL_FALSE;
}

void Adaptivity::disableMovementConnectivity(){
  if(verbose)
    cout<<"void Adaptivity::disableMovementConnectivity()\n";
  TOGTHR = LOGICAL_FALSE;
}

void Adaptivity::disableRefinement(){
  if(verbose)
    cout<<"void Adaptivity::disableRefinement()\n";
  TWOSTG = LOGICAL_TRUE;

  AdaptOpts[0] = LOGICAL_FALSE;
  AdaptOpts[2] = LOGICAL_FALSE;
  AdaptOpts[3] = LOGICAL_FALSE;
  AdaptOpts[4] = LOGICAL_FALSE;
  AdaptOpts[5] = LOGICAL_FALSE;
}

void Adaptivity::disableSurfaceLock(){
  if(verbose)
    cout<<"void Adaptivity::disableSurfaceLock()\n";
  SRFGMY = LOGICAL_FALSE;
}

#ifdef HAVE_VTK
/// Get a vtkUnstructuredGrid object that contains the adapted mesh.
vtkUnstructuredGrid* Adaptivity::get_adapted_vtu(){
  if(verbose)
    cout<<"void Adaptivity::get_adapted_vtk()\n";

  vtkUnstructuredGrid *aug = vtkUnstructuredGrid::New();

  // Set up points.
  vtkPoints *points = vtkPoints::New();
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(newNNodes);
  for(int i=0;i<newNNodes;i++){
    points->SetPoint(i, newX[i], newY[i], newZ[i]);
  }
  aug->SetPoints(points);
  aug->Update();
  points->Delete();
  
  // Set up elements.
  int *newENList = &(intBuffer[NWENLS-1]);
  for(int i=0;i<newNElements;i++){
    vtkIdType pts[4];
    //for(size_t j=0;i<4;j++)
    for(int j=3;j>=0;j--)
      pts[j] = newENList[i*4+j] - 1;
    aug->InsertNextCell(VTK_TETRA, 4, pts);
  }

  // Copy back metric
  vtkDoubleArray *m = vtkDoubleArray::New();
  afloat_t *newMetric = &(floatBuffer[NEWMTX-1]);
  m->SetNumberOfComponents(9);
  m->SetNumberOfTuples(newNNodes);
  m->SetName("metric");
  for(int i=0;i<newNNodes;i++){
    m->SetTuple9(i,
                 newMetric[i*9  ], newMetric[i*9+1], newMetric[i*9+2],
                 newMetric[i*9+3], newMetric[i*9+4], newMetric[i*9+5],
                 newMetric[i*9+6], newMetric[i*9+7], newMetric[i*9+8]);
  }
  aug->GetPointData()->AddArray(m);
  aug->Update();
  m->Delete();

  if(interpolate){
    // Copy interpolated fields.
    afloat_t *new_fields = &(floatBuffer[NEWFLD-1]);
    size_t pos=0;
    for(int n=0;n<nfields;n++){
      vtkDataArray *array = ug->GetPointData()->GetArray(n);
      if(string("metric")==string(array->GetName()))
	continue;
      
      size_t ncomp = array->GetNumberOfComponents();
      
      vtkDoubleArray *m = vtkDoubleArray::New();
      m->SetNumberOfComponents(ncomp);
      m->SetNumberOfTuples(newNNodes);
      m->SetName(array->GetName());

      vector<double> tuple(ncomp);
      for(int i=0;i<newNNodes;i++){
	for(size_t j=0;j<ncomp;j++)
	  tuple[j] = new_fields[pos++];
	m->SetTuple(i, &(tuple[0]));
      }
      aug->GetPointData()->AddArray(m);
      aug->Update();
      m->Delete();
    }
  }
  
  return aug;
}
#endif

void Adaptivity::getHalo(int *_numPrivateNodesint, 
                         int *_ATOSEN, int *_ATOREC, 
                         int *_Gather, int *_Scatter){
  if(verbose)
    cout<<"void Adaptivity::getHalo(int *, int *, int *, int *, int *)\n";

  if(NProcs>1){
    *_numPrivateNodesint = newNPrivateNodes;
    memcpy(_ATOSEN,  &(newATOSEN[0]),  sizeof(int)*(NProcs+1));
    memcpy(_ATOREC,  &(newATOREC[0]),  sizeof(int)*(NProcs+1));
    memcpy(_Gather,  &(newGather[0]),  sizeof(int)*NGather);
    memcpy(_Scatter, &(newScatter[0]), sizeof(int)*NScatter);
  }
}

int Adaptivity::getNProcessors() const{
  if(verbose)
    cout<<"int Adaptivity::getNProcessors() const\n";
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    return MPI::COMM_WORLD.Get_size();
  }
#endif
  return 1;
}

void Adaptivity::getMeshDimensions(int *_NNodes, int *_NElements, int *_NSElements){
  if(verbose)
    cout<<"void Adaptivity::getMeshDimensions(int *, int *, int *)\n";
  *_NNodes     = newNNodes;
  *_NElements  = newNElements;
  *_NSElements = newNSElements;
}

void Adaptivity::get_surface_ids(vector<int> &sids){
  if(verbose)
    cout<<"void Adaptivity::get_surface_ids(vector<int> &sids)\n";
  
  sids.resize(newNSElements);
  for(int i=0;i<newNSElements;i++)
    sids[i] = intBuffer[NWSFID-1+i];
}

void Adaptivity::get_surface_mesh(vector<int> &SENList){
  if(verbose)
    cout<<"void Adaptivity::getSurfaceMesh(vector<int> &SENList)\n";
  
  SENList.resize(3*newNSElements);
  for(int i=0;i<newNSElements*3;i++)
    SENList[i] = intBuffer[NWSNLS-1+i]-1;
}

void Adaptivity::getVertices(afloat_t *_X, afloat_t *_Y, afloat_t *_Z){
  if(verbose)
    cout<<"void Adaptivity::getVertices(afloat_t *, afloat_t *, afloat_t *)\n";
  memcpy(_X, &(floatBuffer[NWNODX-1]), newNNodes*sizeof(afloat_t));
  memcpy(_Y, &(floatBuffer[NWNODY-1]), newNNodes*sizeof(afloat_t));
  memcpy(_Z, &(floatBuffer[NWNODZ-1]), newNNodes*sizeof(afloat_t));
}

void Adaptivity::getVolumeID(int *_newvolumeID){
  if(verbose)
    cout<<"void Adaptivity::getVolumeID(int *)\n";
  memcpy(_newvolumeID, &(intBuffer[NWELRG-1]), newNElements*sizeof(int));
}

void Adaptivity::getVolumeMesh(int *_newENList){
  if(verbose)
    cout<<"void Adaptivity::getVolumeMesh(int *)\n";
  memcpy(_newENList, &(intBuffer[NWENLS-1]), 4*newNElements*sizeof(int));
}

void Adaptivity::setFunctionalTolerance(afloat_t tol){
  if(verbose)
    cout<<"void Adaptivity::setFunctionalTolerance("<<tol<<")\n";
  MESTP1 = tol;
}

void Adaptivity::setHalo(int _numPrivateNodesint, 
                         int *_ATOSEN, int *_ATOREC, 
                         int *_Gather, int *_Scatter){
  if(verbose)
    cout<<"void Adaptivity::setHalo(int, int *, int *, int *, int *)\n";
  NPrivateNodes = _numPrivateNodesint;
  ATOSEN  = _ATOSEN;
  ATOREC  = _ATOREC;
  Gather  = _Gather;
  Scatter = _Scatter;
}

void Adaptivity::set_adapt_sweeps(int _iterations){
  if(verbose)
    cout<<"void Adaptivity::set_adapt_sweeps("<<_iterations<<")\n";
  MaxNumberAdaptIterations = _iterations;
}

void Adaptivity::set_metric(afloat_t *_metric){
  if(verbose)
    cout<<"void Adaptivity::setMetric(afloat_t *)\n";
  size_t NNodes = X.size();
  Metric.resize(9*NNodes);
  for(size_t i=0;i<9*NNodes;i++){
    Metric[i] = _metric[i];
  }
}

void Adaptivity::set_surface_ids(vector<int> &_sids){
  if(verbose)
    cout<<"void Adaptivity::set_surface_ids(vector<int> &_sids)\n";
  
  if(SENList.empty()){
    cerr<<"WARNING: trying to set surface id's but no surface has been specified\n";
    return;
  }
  
  surfID = _sids;
  
  return;
}

void Adaptivity::set_surface_mesh(vector<int> &_SENList){
  if(verbose)
    cout<<"void Adaptivity::set_surface_mesh(vector<int> &_SENList)\n";
  
  SENList = _SENList;
  for(vector<int>::iterator it=SENList.begin();it!=SENList.end();++it)
    *it = *it + 1;

  return;
}

void Adaptivity::set_points(afloat_t *_X, afloat_t *_Y, afloat_t *_Z, int _NNodes){
  if(verbose)
    cout<<"void Adaptivity::set_points(afloat_t *, afloat_t *, afloat_t *, int)\n";
  NPrivateNodes = _NNodes;
  X.resize(_NNodes);
  Y.resize(_NNodes);
  Z.resize(_NNodes);
  for(int i=0;i<_NNodes;i++){
    X[i]=_X[i];
    Y[i]=_Y[i];
    Z[i]=_Z[i];
  }

  return;
}

void Adaptivity::setVolumeID(int *_volumeID){
  if(verbose)
    cout<<"void Adaptivity::setVolumeID(int *)\n";
  volumeID = _volumeID;
}

void Adaptivity::set_volume_mesh(int *_ENList, int _NElements){
  if(verbose)
    cout<<"void Adaptivity::setVolumeMesh(int *, int)\n";

  ENList.resize(4*_NElements);
  for(int i=0;i<4*_NElements;i++)
    ENList[i] = _ENList[i];

  return;
}

#ifdef HAVE_VTK
/// Set the input from a vtkUnstructuredGrid
void Adaptivity::set_from_vtk(vtkUnstructuredGrid *_ug, bool _interpolate){
  if(verbose)
    cout<<"void Adaptivity::set_from_vtk(const vtkUnstructuredGrid *ug, bool interpolate)\n";
  
  ug = _ug;
  interpolate = _interpolate;

  size_t NNodes = ug->GetNumberOfPoints();
  NPrivateNodes = NNodes;
  X.resize(NNodes);
  Y.resize(NNodes);
  Z.resize(NNodes);
  double r[3];
  for(size_t i=0;i<NNodes;i++){
    ug->GetPoints()->GetPoint(i, r);
    X[i] = r[0];
    Y[i] = r[1];
    Z[i] = r[2];
  }
  
  int NElements = ug->GetNumberOfCells();
  ENList.resize(4*NElements);
  for(int i=0;i<NElements;i++){
    assert(ug->GetCell(i)->GetCellType()==VTK_TETRA);
    vtkTetra *tetra = (vtkTetra *)ug->GetCell(i);
    
    //for(size_t j=0;j<4;j++){
    for(int j=3;j>=0;j--){
      ENList[i*4+j] = tetra->GetPointId(j) + 1;
    }
  }

  vtkDataArray *m = ug->GetPointData()->GetArray("metric");
  if(m!=NULL){
    Metric.resize(9*NNodes);
#ifdef DOUBLEP
    for(size_t i=0;i<NNodes;i++){
      m->GetTuple(i, &(Metric[i*9]));      
    }
#else
    double lMetric[9];
    for(size_t i=0;i<NNodes;i++){
      m->GetTuple(i, lMetric);
      for(size_t j=0;j<9;j++){
        Metric[9*i+j]=lMetric[j];
      }
    }
#endif
  }

  if(interpolate){
    nfields = ug->GetPointData()->GetNumberOfArrays();
    for(int i=0;i<ug->GetPointData()->GetNumberOfArrays();i++){
      vtkDataArray *array = ug->GetPointData()->GetArray(i);
      if(string("metric")==string(array->GetName())){
        nfields--;
	continue;
      }

      size_t ncomp = array->GetNumberOfComponents();
      size_t ntuples = array->GetNumberOfTuples();
      nfreedom.push_back(ncomp);
      for(size_t t=0;t<ntuples;t++){
	double tuple[9];
	array->GetTuple(t, tuple);
	for(size_t j=0;j<ncomp;j++)
	  fields.push_back(tuple[j]);
      }
    }
  }

  return;
}
#endif

double Adaptivity::volume(int n1, int n2, int n3, int n4) const{
  if(verbose)
    cout<<"double Adaptivity::volume(int, int, int, int) const\n";

  n1 = n1 - numberingOffset;
  n2 = n2 - numberingOffset;
  n3 = n3 - numberingOffset;
  n4 = n4 - numberingOffset;
  
  return volume(X[n1], Y[n1], Z[n1],
		X[n2], Y[n2], Z[n2],
		X[n3], Y[n3], Z[n3],
		X[n4], Y[n4], Z[n4]);
}

double Adaptivity::volume(int n1, int n2, int n3, afloat_t x, afloat_t y, afloat_t z) const{
  if(verbose)
    cout<<"double Adaptivity::volume(int, int, int, afloat_t, afloat_t, R) const\n";
  n1 = n1 - numberingOffset;
  n2 = n2 - numberingOffset;
  n3 = n3 - numberingOffset;
  
  return volume(X[n1], Y[n1], Z[n1],
		X[n2], Y[n2], Z[n2],
		X[n3], Y[n3], Z[n3],
		x, y, z);
}


double Adaptivity::volume(afloat_t x1, afloat_t y1, afloat_t z1,
                          afloat_t x2, afloat_t y2, afloat_t z2,
                          afloat_t x3, afloat_t y3, afloat_t z3,
                          afloat_t x4, afloat_t y4, afloat_t z4) const{ 
  if(verbose)
    cout<<"double Adaptivity::volume(afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t, afloat_t,y R) const\n";

  double vx1 = x2 - x1;
  double vy1 = y2 - y1;
  double vz1 = z2 - z1;
  
  double vx2 = x3 - x1;
  double vy2 = y3 - y1;
  double vz2 = z3 - z1;
  
  double vx3 = x4 - x1;
  double vy3 = y4 - y1;
  double vz3 = z4 - z1;
  
  double det = 
    vx1*(vy2*vz3 - vy3*vz2) - 
    vx2*(vy1*vz3 - vy3*vz1) + 
    vx3*(vy1*vz2 - vy2*vz1);

  return det/6.0;
}

void Adaptivity::verbose_off(){
  verbose = false;
}

void Adaptivity::verbose_on(){
  verbose = true;
}
