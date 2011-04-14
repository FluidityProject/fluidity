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

#include "DiscreteGeometryConstraints.h"

using namespace std;

bool DiscreteGeometryConstraints::verbose=false;


DiscreteGeometryConstraints::DiscreteGeometryConstraints(){
  if(verbose)
    cout<<"DiscreteGeometryConstraints::DiscreteGeometryConstraints()\n";
  set_coplanar_tolerance(0.9999999);
  ug = NULL;
}

void DiscreteGeometryConstraints::set_coplanar_tolerance(double tol){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::set_coplanar_tolerance("<<tol<<")\n";
  COPLANAR_MAGIC_NUMBER = tol;
}

void DiscreteGeometryConstraints::verbose_off(){
  verbose = false;
}

void DiscreteGeometryConstraints::verbose_on(){
  verbose = true;
}

void DiscreteGeometryConstraints::get_constraints(vector<double> &max_len){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::get_constraints(vector<double> &max_len)\n";

  if(gconstraint.empty())
    find_geometry_constraints();

  size_t NNodes = ug->GetNumberOfPoints();
  max_len.resize(9*NNodes);
  for(size_t i=0; i<NNodes; i++){
    max_len[i*9  ] = gconstraint[i]; max_len[i*9+1] = 0.0;            max_len[i*9+2] = 0.0;
    max_len[i*9+3] = 0.0;            max_len[i*9+4] = gconstraint[i]; max_len[i*9+5] = 0.0;
    max_len[i*9+6] = 0.0;            max_len[i*9+7] = 0.0;            max_len[i*9+8] = gconstraint[i];
  }

  return;
}

void DiscreteGeometryConstraints::get_surface(std::vector<int> &ENList){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::get_surface(vector<int> &ENList)\n";
  ENList = SENList;

  return;
}

void DiscreteGeometryConstraints::set_volume_input(vtkUnstructuredGrid *_ug){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::set_volume_input(vtkUnstructuredGrid *ug)\n";
  ug = _ug;
  ug->Update();

  // Get node-element adjancy list
  size_t NElements = ug->GetNumberOfCells();
  size_t NNodes = ug->GetNumberOfPoints();
  std::deque< std::set<size_t> > NEList(NNodes);
  for(size_t i=0;i<NElements;i++){
    vtkTetra *tetra = (vtkTetra *)ug->GetCell(i);
    for(size_t j=0;j<4;j++){
      NEList[tetra->GetPointId(j)].insert(i);
    }
  }

  // Search for outter surface
  for(size_t i=0;i<NElements;i++){
    vtkTetra *tetra = (vtkTetra *)ug->GetCell(i);
    for(size_t j=0;j<4;j++){
      bool found=false;
      size_t nid1, nid2, nid3;
      if(j==0){
        nid1=tetra->GetPointId(1);
        nid2=tetra->GetPointId(3);
        nid3=tetra->GetPointId(2);
      }else if(j==1){
        nid1=tetra->GetPointId(0);
        nid2=tetra->GetPointId(2);
        nid3=tetra->GetPointId(3);
      }else if(j==2){
        nid1=tetra->GetPointId(0);
        nid2=tetra->GetPointId(3);
        nid3=tetra->GetPointId(1);
      }else{
        nid1=tetra->GetPointId(0);
        nid2=tetra->GetPointId(1);
        nid3=tetra->GetPointId(2);
      }
      
      for(set<size_t>::iterator it=NEList[nid1].begin();it!=NEList[nid1].end();++it){
        if((*it)==i){
          continue;
        }
        
        if(NEList[nid2].find(*it)!=NEList[nid2].end()){
          if(NEList[nid3].find(*it)!=NEList[nid3].end()){
            found=true;
            break;
          }
        }
      }
      if(!found){
        SENList.push_back(nid1);
        SENList.push_back(nid2);
        SENList.push_back(nid3);
      }
    }
  }

  get_coplanar_ids();

  return;
}

void DiscreteGeometryConstraints::set_surface_input(vtkUnstructuredGrid *_ug, vector<int> &_SENList, vector<int> &ids){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::set_surface_input(vtkUnstructuredGrid *ug, vector<int> &_SENList, vector<int> &ids)\n";
  ug = _ug;
  SENList = _SENList;
  
  // Create EEList for surface
  size_t NSElements = SENList.size()/3;
  for(size_t i=0;i<NSElements;i++){
    for(size_t j=0;j<3;j++){
      SNEList[SENList[3*i+j]].insert(i);
    }
  }
  
  coplanar_ids = ids;

  return;
}

#ifdef HAVE_VTK
void DiscreteGeometryConstraints::write_vtk(std::string filename){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::write_vtk(std::string filename)\n";
  
  vtkUnstructuredGrid *sug = vtkUnstructuredGrid::New();
  sug->SetPoints(ug->GetPoints());

  size_t NSElements = SENList.size()/3;
  for(size_t i=0;i<NSElements;i++){
    vtkIdType pts[3];
    for(size_t j=0;j<3;j++){
      pts[j] = SENList[i*3+j];
    }
    sug->InsertNextCell(VTK_TRIANGLE, 3, pts);
  }

  if(coplanar_ids.size()){
    assert(coplanar_ids.size()==NSElements);
    vtkIntArray *field = vtkIntArray::New();
    field->SetNumberOfComponents(1);
    field->SetNumberOfTuples(NSElements);
    field->SetName("coplanar_ids");
    for(size_t i=0;i<NSElements;i++){
      field->SetTuple1(i, coplanar_ids[i]);
    }
    sug->GetCellData()->AddArray(field);
    sug->Update();
    field->Delete();
  }

  size_t NNodes = sug->GetNumberOfPoints();
  if(gconstraint.size()){
    vtkDoubleArray *field = vtkDoubleArray::New();
    field->SetNumberOfComponents(1);
    field->SetNumberOfTuples(NNodes);
    field->SetName("max_lengths");
    for(size_t i=0;i<NNodes;i++){
      field->SetTuple1(i, gconstraint[i]);
    }
    sug->GetPointData()->AddArray(field);
    sug->Update();
    field->Delete();
  }

  vtkXMLUnstructuredGridWriter *sug_writer = vtkXMLUnstructuredGridWriter::New();
  sug_writer->SetFileName(filename.c_str());
  sug_writer->SetInput(sug);
  sug_writer->Write();
  sug_writer->Delete();
  
  sug->Delete();
  return;
}
#endif

void DiscreteGeometryConstraints::normal(int v1, int v2, int v3, double *n){
  double r1[3], r2[3], r3[3];
  ug->GetPoint(v1, r1);
  ug->GetPoint(v2, r2);
  ug->GetPoint(v3, r3);
    
  double ax = r3[0] - r2[0];
  double ay = r3[1] - r2[1];
  double az = r3[2] - r2[2];

  double bx = r1[0] - r2[0];
  double by = r1[1] - r2[1];
  double bz = r1[2] - r2[2];
  
  n[0] = (ay * bz - az * by);
  n[1] = (az * bx - ax * bz);
  n[2] = (ax * by - ay * bx);

  double length = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  
  n[0]/=length;
  n[1]/=length;
  n[2]/=length;
  
  return;
}

void DiscreteGeometryConstraints::get_coplanar_ids(){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::get_coplanar_ids()\n";

  // Calculate all element normals
  size_t NSElements = SENList.size()/3;
  vector<double> normals(NSElements*3);
  for(size_t i=0;i<NSElements;i++){
    normal(SENList[3*i], SENList[3*i+1], SENList[3*i+2], &(normals[i*3]));
  }

  // Create EEList for surface
  for(size_t i=0;i<NSElements;i++){
    for(size_t j=0;j<3;j++){
      SNEList[SENList[3*i+j]].insert(i);
    }
  }
  
  EEList.resize(NSElements*3);
  for(size_t i=0;i<NSElements;i++){
    for(size_t j=0;j<3;j++){
      size_t nid1=SENList[i*3+(j+1)%3];
      size_t nid2=SENList[i*3+(j+2)%3];
      for(set<size_t>::iterator it=SNEList[nid1].begin();it!=SNEList[nid1].end();++it){
        if((*it)==i){
          continue;
        }
        
        if(SNEList[nid2].find(*it)!=SNEList[nid2].end()){
          EEList[i*3+j] = *it;
          break;
        }
      }
    }
  }

  // Form patches
  coplanar_ids.resize(NSElements);
  for(vector<int>::iterator it=coplanar_ids.begin(); it!=coplanar_ids.end(); ++it)
    *it = 0;
  
  size_t current_id = 1;
  for(size_t pos = 0;pos<NSElements;){
    // Create a new starting point
    double *ref_normal=NULL;
    for(size_t i=pos;i<NSElements;i++){{
        if(coplanar_ids[i]==0){
          // This is the first element in the new patch
          pos = i;
          coplanar_ids[pos] = current_id;
          ref_normal = &(normals[pos*3]);
          break;
        }
      }
      
      // Jump out of this while loop if we are finished
      if(i==NSElements)
        break;
    }
    
    // Initialise the front
    set<int> front;
    front.insert(pos);
          
    // Advance this front
    while(!front.empty()){
      int sele = *front.begin();
      front.erase(front.begin());

      // Check surrounding surface elements:      
      for(int i=0; i<3; i++){
        int sele2 = EEList[sele*3+i];
        
        if(coplanar_ids[sele2]>0)
          continue;
          
        double coplanar =
          ref_normal[0]*normals[sele2*3  ]+
          ref_normal[1]*normals[sele2*3+1]+
          ref_normal[2]*normals[sele2*3+2];

        if(coplanar>=COPLANAR_MAGIC_NUMBER){
          front.insert(sele2);
          coplanar_ids[sele2] = current_id;
        }
      }
    }
    current_id++;
    pos++;
  }

  return;
}

void DiscreteGeometryConstraints::get_coplanar_ids(std::vector<int> &ids){
  if(verbose)
    cout<<"void get_coplanar_ids(std::vector<int> &ids)\n";

  if(coplanar_ids.empty())
    get_coplanar_ids();

  ids = coplanar_ids;
}

void DiscreteGeometryConstraints::find_geometry_constraints(){
  if(verbose)
    cout<<"void DiscreteGeometryConstraints::find_geometry_constraints()\n";

  // Find the maximum lengths that will be used for constraining
  // the metric.
  double bbox[6];
  ug->GetBounds(bbox);
  double max_dist = max(max(bbox[1]-bbox[0], bbox[3]-bbox[2]), bbox[5]-bbox[4]);
    
  // Initialise geometry constraints
  size_t NNodes = ug->GetNumberOfPoints();
  gconstraint.resize(NNodes);
  for(size_t i=0;i<NNodes;i++){
    gconstraint[i] = max_dist;
  }

  size_t NSElements=coplanar_ids.size();
  size_t npatches=0;
  for(size_t i=0;i<NSElements;i++){
    npatches = max((int)npatches, coplanar_ids[i]);
  }

  // Loop over each patch, merging in the constraints from each patch
  // into the nodes that lie along geometry edges.
  for(size_t p=0;p<npatches;p++){
    
    map<size_t, set<size_t> > corner;
    // Identify the corner nodes of this patch.
    for(size_t ele=0;ele<NSElements;ele++){
      if(coplanar_ids[ele]!=(int)p)
        continue;
      
      for(size_t j=0;j<3;j++){
        size_t nid = SENList[ele*3+j];
        
        set<size_t> patches;
        for(set<size_t>::iterator it=SNEList[nid].begin();it!=SNEList[nid].end();++it){
          patches.insert(coplanar_ids[*it]);
        }
        
        if(patches.size()>2){
          corner[nid] = patches;
        }
      }
    }

    // I'm taking a short cut from the full algorithm - it will serve
    // as a good approximation when the surface is complex.
    for(map<size_t, set<size_t> >::iterator it=corner.begin(); it!=corner.end();++it){
      double dist = max_dist;
      int n0 = it->first;
      for(map<size_t, set<size_t> >::iterator jt=corner.begin(); jt!=corner.end();++jt){
        int n1 = jt->first;
        if(n0==n1)
          continue;

	double r0[3], r1[3];
	ug->GetPoint(n0, r0);
	ug->GetPoint(n1, r1);

        double dx = r0[0] - r1[0];
        double dy = r0[1] - r1[1];
        double dz = r0[2] - r1[2];
        
        double d = sqrt(dx*dx+dy*dy+dz*dz);

        dist = min(dist, d);
      }
      gconstraint[n0] = min(dist, gconstraint[n0]);
    }
  }

  return;
}
