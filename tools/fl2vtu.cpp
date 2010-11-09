/*  Copyright (C) 2006 Imperial College London and others.

    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk

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

#include <confdefs.h>

#include <iostream>
#ifdef HAVE_VTK

#include <algorithm>
#include <cassert>
#include <fstream>
#include <deque>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <set>

#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
}

using namespace std;

bool verbose, legacy, periodic;

void usage(char *name){
  cout<<name<<" [-l] [-z] <fluidity-project-name> <dump-id> [<last-dump-id>]\n";
  cout<<"\nConvert a fluidity dump file (or range of) into a VTK unstructured grid XML file(s). "
      <<"Serial fluidity dumps are named <fluidity-project-name>_<dump-id>.vtu "
      <<"while parallel fluidity dumps consist of a metafile named "
      <<"<fluidity-project-name>_<dump-id>.pvtu, that has all the partition "
      <<"information, and an individual .vtu file for each of the partitions "
      <<"in the domains, <fluidity-project-name>_<dump-id>_<partition number>.vtu\n"
      <<"\n-h prints this\n"
      <<"\n-l is required for old parallel dump files\n"
      <<"\n-m merge dataset into a single file if parallel\n"
      <<"\n-v verbose\n";
  return;
}

// Split a line up into tokens
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " "){
  tokens.clear();
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos){
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

bool logical2bool(char flogical){
  return (flogical=='t')||(flogical=='T');
}


// Convert word to uppercase.
void uppercase(string& word){
  for(string::iterator c = word.begin(); c != word.end(); ++c){
    if( (*c >= 'a')&&(*c <= 'z') ){
      *c = (*c - 'a') + 'A';
    }
  }
}

// Read array from fluidity dump file
int ReadArray(ifstream &file, float *array, int cnt){
  for(int i=0;i<cnt;i++){
    file>>array[i];
    assert(!file.eof());
  }
  return(0);
}

// Read array from fluidity dump file
int ReadArray(ifstream &file, int *array, int cnt){
  for(int i=0;i<cnt;i++){
    file>>array[i];
    assert(!file.eof());
  }
  return(0);
}

// Get dump-file names
int get_dump_file_names(string project_name, deque<string> &dump_file_names, int dump_id){
  char dump_name[2048];
  sprintf(dump_name,"%s.d.%d",project_name.c_str(),dump_id);

  ifstream file(dump_name, ios::in);
  if(file.is_open()){
    dump_file_names.push_back(string(dump_name));
    file.close();
  }else{
    // So...must be parallel dump
    for(size_t partition=1;;partition++){
      sprintf(dump_name,"%s.d%d.%d",project_name.c_str(),partition,dump_id);
      file.open(dump_name, ios::in);
      if(file.is_open()){
        dump_file_names.push_back(string(dump_name));
        file.close();
      }else{
        break;
      }
    }
  }

  return dump_file_names.size();
}

int ReadDump(string dump_name,
             size_t &n_private_nodes, size_t &nelements, int &vtk_cell_type,
             vector<float> &X, vector<float> &Y, vector<float> &Z,
             vector<int> &NDGLNO,
             vector< vector<float> > &Ux, vector< vector<float> > &Uy, vector< vector<float> > &Uz,
             vector< vector<float> > &data_arrays,
             vector<float> &Pressure,
             map<unsigned, unsigned> &periodic_renumbering,
             bool merge, vector<int> &COLGA, vector<int> &SCATER, vector<int> &ATOSEN, vector<int> &ATOREC){
  
  ifstream file(dump_name.c_str(), ios::in);
  assert(file.is_open());
  
  char buffer[2048];
  vector<string> tokens;
  
  // Read header information
  int SAVRES, NTSOL, NLOC, NPHASE;
  file>>SAVRES>>NTSOL>>NLOC>>NPHASE;
  if(verbose) cout<<"SAVRES, NTSOL, NLOC, NPHASE = "<<SAVRES<<" "<<NTSOL<<" "<<NLOC<<" "<<NPHASE<<endl;
  
  int nnodes, NNODPP;
  file>>nnodes>>nelements>>n_private_nodes>>NNODPP;
  if(verbose) cout<<"nnodes, nelements, n_private_nodes, NNODPP = "<<nnodes<<" "<<nelements<<" "<<n_private_nodes<<" "<<NNODPP<<endl;
  
  int data_length = nnodes;
  if(legacy)
    data_length = n_private_nodes;
  
  int FREDOP;
  file>>FREDOP;
  
  int data_length2 = FREDOP;
  if(legacy)
    data_length2 = NNODPP;
  
  int MLOC = NLOC;
  if(SAVRES>=100)
    file>>MLOC;
  MLOC=abs(MLOC);
  
  int VERSIO,XNONOD;
  if(NLOC<0){
    file>>VERSIO>>XNONOD;
    
    // read return character
    file.getline(buffer, 2048);

    file.getline(buffer, 2048);

    Tokenize(buffer, tokens, " ");
    // assert(tokens.size()==1); - this line contains addition logicals which can be ignored
    periodic=logical2bool(tokens[0][0]);
    NLOC = abs(NLOC);
  }else{
    VERSIO=0;
    XNONOD=nnodes;
    periodic=false;
  }
  
  if(periodic)
    cerr<<"WARNING: No support written for periodic boundaries\n";

  if(verbose) cout<<"VERSIO, XNONOD, periodic = "<<VERSIO<<" "<<XNONOD<<" "<<periodic<<endl;

  float ACCTIM, DT, D0, R0;
  file>>ACCTIM>>DT>>D0>>R0;
  if(verbose) cout<<"ACCTIM, DT, D0, R0 = "<<ACCTIM<<" "<<DT<<" "<<D0<<" "<<R0<<endl;

  // read return character
  file.getline(buffer, 2048);

  bool NAV, D3;
  file.getline(buffer, 2048);
  Tokenize(buffer, tokens, " ");
  assert(tokens.size()==2);
  
  NAV=logical2bool(tokens[0][0]);
  D3=logical2bool(tokens[1][0]);
  if(verbose) cout<<"NAV,D3 = "<<NAV<<" "<<D3<<endl;
  
  vector<int> TELEDI(abs(NTSOL), 0);
  if(NTSOL<0){
    if(verbose) cout<<"reading TELEDI[]\n";
    NTSOL=abs(NTSOL);
    ReadArray(file, &(TELEDI[0]), NTSOL);
  }

  // get arrays
  if(verbose) cout<<"getting data arrays\n";
  // int nscalars = max((int)NTSOL, (int)Scalars.size());
  int nscalars = abs(NTSOL);
  data_arrays.resize(nscalars);
  vector< vector<float> > ETNEW(nscalars);
  for(int i=0;i<nscalars;i++){
    data_arrays[i].resize(max(nnodes,XNONOD));
    ReadArray(file, &(data_arrays[i][0]), data_length);
    if(TELEDI[i]==1){
      if(verbose) cout<<"reading ETNEW\n";
      ETNEW[i].resize(max(nnodes,XNONOD));
      ReadArray(file, &(ETNEW[i][0]), data_length2);
    }
  }

  if(NAV){
    Ux.resize(max(1,NPHASE));
    Uy.resize(max(1,NPHASE));
    Uz.resize(max(1,NPHASE));
    for(int i=0;i<max(1,NPHASE);i++){
      Ux[i].resize(nnodes);
      ReadArray(file, &(Ux[i][0]), data_length);
      
      Uy[i].resize(nnodes);
      ReadArray(file, &(Uy[i][0]), data_length);
      
      Uz[i].resize(nnodes);
      if(D3){
        ReadArray(file, &(Uz[i][0]), data_length);
      }else{
        for(int j=0;j<nnodes;j++)
          Uz[i][j] = 0;
      }
    }
    
    if(SAVRES>=2){
      if(verbose) cout<<"Reading Pressure\n";
      Pressure.resize(FREDOP);
      ReadArray(file, &(Pressure[0]), data_length2);
    }
    
    vector<float> FRAPHA;
    if(NPHASE>1){
      if(verbose) cout<<"Reading FRAPHA\n";
      FRAPHA.resize(FREDOP*NPHASE);
      ReadArray(file, &(FRAPHA[0]), FREDOP*NPHASE);
    }
  }

  bool YES = false;
  for(int i=0;i<NTSOL;i++)
    if(TELEDI[i]!=0)
      YES = true;

  vector< vector<float> > EDNEW;
  if(YES&&NAV&&(NPHASE>1)&&(SAVRES>3)){
    EDNEW.resize(max(1,NPHASE));
    for(int i=0;i<max(1,NPHASE);i++){
      EDNEW[i].resize(NNODPP);
      ReadArray(file, &(EDNEW[i][0]), NNODPP);
    }
  }

  NDGLNO.resize(nelements*NLOC);
  vector<int> XONDGL;
  vector<int> PNDGLN(nelements*MLOC);
  if(SAVRES>=3){
    if(verbose) cout<<"getting coordinates\n";
    X.resize(XNONOD);
    ReadArray(file, &(X[0]), XNONOD);

    Y.resize(XNONOD);
    ReadArray(file, &(Y[0]), XNONOD);

    Z.resize(XNONOD);
    if(D3){
      ReadArray(file, &(Z[0]), XNONOD);
    }else{
      for(int i=0;i<XNONOD;i++)
        Z[i] = 0.0;
    }

    if(SAVRES>=4){
      if(verbose) cout<<"reading regular mesh\n";
      ReadArray(file, &(NDGLNO[0]), nelements*NLOC);

      if(periodic){
        if(verbose) cout<<"reading periodic mesh\n";
        XONDGL.resize(nelements*NLOC);
        ReadArray(file, &(XONDGL[0]), nelements*NLOC);
        for(size_t i=0;i<nelements;i++)
          for(int j=0;j<NLOC;j++)
            periodic_renumbering[XONDGL[i*NLOC+j]-1] = NDGLNO[i*NLOC+j]-1;
        NDGLNO.swap(XONDGL);
        XONDGL.clear();
        nnodes = XNONOD;
      }

    }else{
      cerr<<"Don't know what to do when there is no mesh! (SAVRES<4)\n";
      exit(-1);
    }

    // Sort out fortran numbering
    for(vector<int>::iterator it=NDGLNO.begin();it!=NDGLNO.end();it++){
      *it = *it - 1;
      assert(*it>=0);
      assert(*it<nnodes);
    }
  }

  if(NAV&&(SAVRES>=100)){
    if(verbose) cout<<"reading pressure mesh"<<endl;
    ReadArray(file, &(PNDGLN[0]), PNDGLN.size());
  }
  
  // but we don't use
  PNDGLN.clear();
  
  if((NLOC==4)&&D3){
    vtk_cell_type = VTK_TETRA;
  }else if((NLOC==8)&&D3){
    vtk_cell_type = VTK_HEXAHEDRON;
  }else if((NLOC==4)&&(!D3)){
    vtk_cell_type = VTK_QUAD;
  }else{
    cerr<<"unsupported cell type - Gerard went for beer before finishing\n";
    cerr<<NLOC<<" "<<D3<<endl;
    exit(-1);
  }
  
  if(n_private_nodes!=X.size()){
    if(merge){
      int NPROCS=0, NCOLGA=0, NSCATE=0;
      file>>NPROCS;
      if(!file.eof()){
  assert(NPROCS>1);
  file>>NCOLGA;
  assert(!file.eof());
  file>>NSCATE;
  assert(!file.eof());
  
  COLGA.resize(NCOLGA);
  ReadArray(file, &(COLGA[0]), NCOLGA);
  SCATER.resize(NSCATE);
  ReadArray(file, &(SCATER[0]), NSCATE);
  ATOSEN.resize(NPROCS+1);
  ReadArray(file, &(ATOSEN[0]), NPROCS+1);
  ATOREC.resize(NPROCS+1);
  ReadArray(file, &(ATOREC[0]), NPROCS+1);
      }
    }else{
      // Strip halo 2
      set<size_t> required_halo_nodes;
      for(size_t i=0;i<nelements;i++){
  bool is_halo2=true;
  for(size_t j=0;j<NLOC;j++){
    if(NDGLNO[i*NLOC+j]<n_private_nodes){
      is_halo2=false;
      break;
    }
  }
  
  if(!is_halo2){
    for(size_t j=0;j<NLOC;j++){
      if(NDGLNO[i*NLOC+j]>=n_private_nodes){
        required_halo_nodes.insert(NDGLNO[i*NLOC+j]);
      }
    }
  }
      }
      
      map<size_t, size_t> renumbering;
      size_t next_nid=n_private_nodes;
      for(set<size_t>::const_iterator it=required_halo_nodes.begin(); it!=required_halo_nodes.end(); it++){
  renumbering[*it] = next_nid++;
      }
      
      map<size_t, size_t> erenumbering;
      size_t new_nelements=0;
      for(size_t i=0;i<nelements;i++){
  bool keep_element=true;
  for(size_t j=0;j<NLOC;j++){
    size_t nid=NDGLNO[i*NLOC+j];
    if(nid<n_private_nodes){
      NDGLNO[new_nelements*NLOC+j] = nid;
    }else{
      if(renumbering.find(nid)!=renumbering.end()){
        NDGLNO[new_nelements*NLOC+j] = renumbering[nid];
      }else{
        keep_element=false;
        break;
      }
    }
  }
  if(keep_element){
    erenumbering[i] = new_nelements++;
  }
      }
      nelements = new_nelements;
      NDGLNO.resize(nelements*NLOC);

      for(map<size_t, size_t>::const_iterator it=renumbering.begin(); it!=renumbering.end(); it++){
  X[it->second] = X[it->first];
  Y[it->second] = Y[it->first];
  Z[it->second] = Z[it->first];
  
  for(size_t j=0;j<Ux.size();j++)
    Ux[j][it->second] = Ux[j][it->first];

  for(size_t j=0;j<Uy.size();j++)
    Uy[j][it->second] = Uy[j][it->first];
  
  for(size_t j=0;j<Uz.size();j++)
    Uz[j][it->second] = Uz[j][it->first];
  
  for(size_t j=0;j<data_arrays.size();j++)
    data_arrays[j][it->second] = data_arrays[j][it->first];  
      }
      
      if(Pressure.size()==X.size()){
  for(map<size_t, size_t>::const_iterator it=renumbering.begin(); it!=renumbering.end(); it++){
    Pressure[it->second] = Pressure[it->first];  
  }
  Pressure.resize(n_private_nodes+renumbering.size());
      }else{
  for(map<size_t, size_t>::const_iterator it=erenumbering.begin(); it!=erenumbering.end(); it++){
    Pressure[it->second] = Pressure[it->first];
  }
  Pressure.resize(erenumbering.size());
      }
      
      size_t new_size=n_private_nodes+renumbering.size();
      
      X.resize(new_size);
      Y.resize(new_size);
      Z.resize(new_size);
  
      for(size_t j=0;j<Ux.size();j++)
  Ux[j].resize(new_size);
      
      for(size_t j=0;j<Uy.size();j++)
  Uy[j].resize(new_size);
      
      for(size_t j=0;j<Uz.size();j++)
  Uz[j].resize(new_size);
      
      for(size_t j=0;j<data_arrays.size();j++)
  data_arrays[j].resize(new_size);

    }
  }
  
  return 0;
}

int Create_vtkUnstructuredGrid(vtkUnstructuredGrid *dataSet, int partition_number,
                               map<string, int> &Scalars,
                               map<string, deque<int> > &Vectors,
                               map<string, deque<int> > &Tensors,
                               map<string, deque<int> > &Tensors2,
                               size_t n_private_nodes, size_t nelements,   int vtk_cell_type,
                               vector<float> &X, vector<float> &Y, vector<float> &Z,
                               vector<int> &NDGLNO,
                               vector< vector<float> > &Ux, vector< vector<float> > &Uy, vector< vector<float> > &Uz,
                               vector< vector<float> > &data_arrays,
                               vector<float> &Pressure,
                               map<unsigned, unsigned> &periodic_renumbering){

  // Point definitions
  size_t nnodes = X.size();
  if(verbose) cout<<"Point definitions: "<<nnodes<<endl;
  assert(nnodes==Y.size());
  assert(nnodes==Z.size());
  int NLOC = NDGLNO.size()/nelements;
  vtkPoints *newPts = vtkPoints::New();

  for(size_t i=0; i<nnodes; i++)
    newPts->InsertNextPoint(X[i], Y[i], Z[i]);
  dataSet->SetPoints(newPts);
  newPts->Delete();

  if(verbose) cout<<"cell definitions: "<<nelements<<endl;
  vtkIdType Cell[8];

  for(size_t i=0; i<nelements; i++){
    if(vtk_cell_type==VTK_TETRA){
      for(int j=0;j<4;j++)
        Cell[j] = NDGLNO[i*4+j];
      dataSet->InsertNextCell(VTK_TETRA, 4, Cell);
    }else if(vtk_cell_type==VTK_HEXAHEDRON){
      Cell[0] = NDGLNO[i*8+0];
      Cell[1] = NDGLNO[i*8+1];
      Cell[2] = NDGLNO[i*8+3];
      Cell[3] = NDGLNO[i*8+2];
      Cell[4] = NDGLNO[i*8+4];
      Cell[5] = NDGLNO[i*8+5];
      Cell[6] = NDGLNO[i*8+7];
      Cell[7] = NDGLNO[i*8+6];
      dataSet->InsertNextCell(VTK_HEXAHEDRON, 8, Cell);
    }else if(vtk_cell_type==VTK_QUAD){
      Cell[0] = NDGLNO[i*4+0];
      Cell[1] = NDGLNO[i*4+1];
      Cell[2] = NDGLNO[i*4+3];
      Cell[3] = NDGLNO[i*4+2];
      dataSet->InsertNextCell(VTK_QUAD, 4, Cell);
    }else{
      cerr<<"ERROR: Element not recognised: "<<vtk_cell_type<<endl; 
      exit(-1);
    }
  }

  if(verbose) cout<<"scalar definitions\n";
  for(map<string, int>::iterator it=Scalars.begin();it!=Scalars.end();it++){
    string tag;
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(nnodes);
    if((it->first.compare("presur")==0)||(it->first.compare("pressure")==0)){
      if(verbose) cout<<"pressure\n";
      tag = "Pressure";
      if(Pressure.size()==nnodes){
        for(size_t i=0;i<nnodes;i++){
          size_t nid = i;
          if(periodic)
            nid = periodic_renumbering[i];
          newScalars->InsertValue(i, Pressure[nid]);
        }
      }else if(Pressure.size()==nelements){
        newScalars->SetNumberOfTuples(nelements);
        for(size_t i=0;i<nelements;i++)
          newScalars->InsertValue(i, Pressure[i]);
      }else{
        cerr<<"Don't know how to deal with pressure (Pressure.size()="<<Pressure.size()<<", nelements="<<nelements<<", nnodes="<<nnodes<<")\n";
      }
    }else{
      if(it->first.compare("temp")==0){
        if(verbose) cout<<"temperature\n";
        tag = "Temperature";
      }else{
        if(verbose) cout<<it->first<<endl;
        tag = it->first;
      }
      for(size_t i=0;i<nnodes;i++){
        unsigned nid = i;
        if(periodic)
          nid = periodic_renumbering[i];
        newScalars->InsertValue(i, data_arrays[it->second][nid]);
      }
    }
    if(verbose) cout<<"VTKed scalars\n";
    newScalars->SetName(tag.c_str());
    if(((it->first.compare("presur")==0)||(it->first.compare("pressure")==0))&&(Pressure.size()==nelements)){
      dataSet->GetCellData()->AddArray(newScalars);
      dataSet->GetCellData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::SCALARS);
    }else{
      dataSet->GetPointData()->AddArray(newScalars);
      dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::SCALARS);
    }
    newScalars->Delete();
    if(verbose) cout<<"Happy with scalars\n";
  }

  if(verbose) cout<<"cells...\n";
  // Scalar field for domain decomposition
  if(partition_number!=-1){
    vector<int> decomposition(nelements, partition_number);
    for(size_t i=0;i<nelements;i++){
      for(int j=0;j<NLOC;j++){
        if(NDGLNO[i*NLOC+j]>(int)n_private_nodes){
          decomposition[i] = -1;
        }
      }
    }
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(nelements);
    for(size_t i=0;i<nelements;i++)
      newScalars->InsertValue(i, decomposition[i]);
    newScalars->SetName("Partitioning");
    dataSet->GetCellData()->AddArray(newScalars);
    dataSet->GetCellData()->SetActiveAttribute("Partitioning", vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
  }

  if(verbose) cout<<"Velocity vectors"<<endl;
  for(size_t i=0;i<Ux.size();i++){
    string tag;
    vtkFloatArray *newVectors = vtkFloatArray::New();
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(nnodes);
    tag = "Velocity";
    if(Ux.size()>1){
      if(verbose) cout<<"Write phase "<<i+1<<endl;
      char str[20];
      sprintf(str, " Phase %d", i+1);
      tag.append(str);
    }
    for(size_t n=0; n<nnodes; n++){
      int nid = n;
      if(periodic)
        nid = periodic_renumbering[n];
      newVectors->SetTuple3(n,
                            Ux[i][nid],
                            Uy[i][nid],
                            Uz[i][nid]);
    }
    newVectors->SetName(tag.c_str());
    dataSet->GetPointData()->AddArray(newVectors);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::VECTORS);
    newVectors->Delete();
  }

  for(map<string, deque<int> >::iterator it=Vectors.begin();it!=Vectors.end();it++){
    if(it->first.compare(0, 4, "velu")==0)
      continue;

    string tag;
    vtkFloatArray *newVectors = vtkFloatArray::New();
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(nnodes);

    tag = it->first;
    for(size_t n=0; n<nnodes; n++){
      int nid = n;
      if(periodic){
        nid = periodic_renumbering[n];
      }
      if(it->second.size()==2)
        newVectors->SetTuple3(n,
                              data_arrays[it->second[0]][nid],
                              data_arrays[it->second[1]][nid],
                              0.0);
      else
        newVectors->SetTuple3(n,
                              data_arrays[it->second[0]][nid],
                              data_arrays[it->second[1]][nid],
                              data_arrays[it->second[2]][nid]);
    }
    newVectors->SetName(tag.c_str());
    dataSet->GetPointData()->AddArray(newVectors);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::VECTORS);
    newVectors->Delete();
  }

  if(verbose) cout<<"Tensors\n";
  for(map<string, deque<int> >::iterator it=Tensors.begin();it!=Tensors.end();it++){
    vtkFloatArray *newTensors = vtkFloatArray::New();
    newTensors->SetNumberOfComponents(9);
    newTensors->SetNumberOfTuples(nnodes);
    string tag = it->first;
    for(size_t n=0;n<nnodes;n++){
      int nid = n;
      if(periodic){
        nid = periodic_renumbering[n];
      }
      newTensors->SetTuple9(n,
                            data_arrays[it->second[0]][nid], data_arrays[it->second[1]][nid], data_arrays[it->second[2]][nid],
                            data_arrays[it->second[3]][nid], data_arrays[it->second[4]][nid], data_arrays[it->second[5]][nid],
                            data_arrays[it->second[6]][nid], data_arrays[it->second[7]][nid], data_arrays[it->second[8]][nid] );
    }
    newTensors->SetName(tag.c_str());
    dataSet->GetPointData()->AddArray(newTensors);
    dataSet->GetPointData()->SetActiveAttribute(tag.c_str(), vtkDataSetAttributes::TENSORS);
    newTensors->Delete();
  }
  for(map<string, deque<int> >::iterator it=Tensors2.begin();it!=Tensors2.end();it++){
    cerr<<"("<<__FILE__<<", "<<__LINE__<<") - gerard went for beer\n";
  }

  return 0;
}

int main(int argv, char **argc){
  if(argv==1){
    usage(argc[0]);
    exit(-1);
  }

  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;
  char c;
  int getopt_c;
  map<char, string> args;
  while ((getopt_c = getopt(argv, argc, "hlmv")) != EOF){
    c = (char)getopt_c;
    if (c != '?'){
      if (optarg == NULL){
        args[c] = "true";
      }else{
        args[c] = optarg;
      }
    }else{
      if (isprint(optopt)){
        cerr << "Unknown option " << optopt << endl;
      }else{
        cerr << "Unknown option " << hex << optopt << endl;
      }
      usage(argc[0]);
      exit(-1);
    }
  }

  // Help?
  if(args.find('h')!=args.end()){
    usage(argc[0]);
    exit(-1);
  }

  // Check if legacy behaviour is required
  legacy=(args.find('l')!=args.end());

  // Verbose?
  verbose=(args.find('v')!=args.end());

  if(((argv-optind)!=2)&&
     ((argv-optind)!=3)){
    usage(argc[0]);
    exit(-1);
  }

  string project_name;
  if(argc[argv-2][0]=='/'){
    project_name = argc[optind];
  }else{
    project_name = getenv("PWD");
    project_name+="/";
    project_name+=argc[optind];
  }
  if(verbose) cout<<"project name = "<<project_name<<endl;

  // Field data storage
  map<string, int> Scalars;
  map<string, deque<int> > Vectors;
  map<string, deque<int> > Tensors;
  map<string, deque<int> > Tensors2;

  // Read array description file
  {
    string fld_file_name(project_name);
    fld_file_name+=".fld";
    ifstream fld_file(fld_file_name.c_str(), ios::in);
    if(!fld_file.is_open()){
      cerr<<"field description file not found\n";
      exit(-1);
    }

    string buffer;
    vector<string> tokens;
    for(;;){
      if(getline(fld_file, buffer)){
  Tokenize(buffer.c_str(), tokens, " ");
  unsigned tcnt = tokens.size();
  
  if(tcnt<1)
    continue;
  
  if( (tokens[0][0] == '#') || (tokens[0][0] == '@') )
    continue;
  
  if(tcnt<3){
    cerr<<"corrupt field description file\n";
    exit(-1);
  }
  
  if(tokens[0].compare("SCALAR")==0){
    assert(tokens.size()==3);
    Scalars[tokens[1]] = atoi(tokens[2].c_str())-1;
  }else if(tokens[0].compare("VECTOR")==0){
    assert((tokens.size()==4)||(tokens.size()==5));
    for(size_t i=0;i<tokens.size()-2;i++)
      Vectors[tokens[1]].push_back(atoi(tokens[2+i].c_str())-1);
  }else if(tokens[0].compare("TENSOR")==0){
    assert(tokens.size()==11);
    for(int i=0;i<9;i++){
      Tensors[tokens[1]].push_back(atoi(tokens[2+i].c_str())-1);
    }
  }else if(tokens[0].compare("CONDUCTIVITY")==0){
    assert(tokens.size()==8);
    for(int i=0;i<6;i++){
      Tensors2[tokens[1]].push_back(atoi(tokens[2+i].c_str())-1);
    }
  }else{
    cerr<<"corrupt field description file\n";
    exit(-1);
  }
      }else{
        break;
      }
    }
  }

  int first_dump_id = atoi(argc[optind+1]);
  int last_dump_id = first_dump_id;
  if((argv-optind)==3)
    last_dump_id = atoi(argc[optind+2]);
  for(int dump_id = first_dump_id; dump_id<=last_dump_id;dump_id++){

    // Get dump-file names
    deque<string> dump_file_names;
    size_t npartitions = get_dump_file_names(project_name, dump_file_names, dump_id);

    deque<size_t> n_private_nodes(npartitions), nelements(npartitions);
    int vtk_cell_type;
    deque< vector<float> > X(npartitions), Y(npartitions), Z(npartitions);
    deque< vector<int> > NDGLNO(npartitions);
    deque< vector< vector<float> > > Ux(npartitions), Uy(npartitions), Uz(npartitions);
    deque< vector< vector<float> > > data_arrays(npartitions);
    deque< vector<float> > Pressure(npartitions);
    deque< map<unsigned, unsigned> > periodic_renumbering(npartitions);
    deque< vector<int> > COLGA(npartitions), SCATER(npartitions), ATOSEN(npartitions), ATOREC(npartitions);
    for(size_t part=0; part<npartitions; part++){
      ReadDump(dump_file_names[part],
               n_private_nodes[part], nelements[part], vtk_cell_type,
               X[part], Y[part], Z[part],
               NDGLNO[part],
               Ux[part], Uy[part], Uz[part],
               data_arrays[part], Pressure[part],
               periodic_renumbering[part],
               args.find('m')!=args.end(), COLGA[part], SCATER[part], ATOSEN[part], ATOREC[part]);

      // We must get global numbering information first if we're merging datasets.
      // deque< deque<int> > global_numbering(npartitions);

      // get halo and numbering information
      // get it to compile at this point before going on.

      if((args.find('m')==args.end())&&(npartitions>1)){
        if(verbose) cout<<"writting parallel file\n";

        // Construct vtk object
        if(verbose) cout<<"Construct vtk object\n";
        vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
        Create_vtkUnstructuredGrid(dataSet, part,
                                   Scalars, Vectors, Tensors, Tensors2,
                                   n_private_nodes[part], nelements[part], vtk_cell_type,
                                   X[part], Y[part], Z[part],
                                   NDGLNO[part],
                                   Ux[part], Uy[part], Uz[part],
                                   data_arrays[part], Pressure[part],
                                   periodic_renumbering[part]);

        char vtk_file_name[2048];
        sprintf(vtk_file_name,"%s_%d.pvtu",project_name.c_str(),dump_id);
        vtkXMLPUnstructuredGridWriter *writer= vtkXMLPUnstructuredGridWriter::New();
        vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
        writer->SetDataModeToBinary();
        writer->SetFileName(vtk_file_name);
        writer->SetNumberOfPieces(npartitions);
        writer->SetGhostLevel(1);
        writer->SetStartPiece(part);
        writer->SetEndPiece(part);
        writer->SetInput(dataSet);
        writer->SetCompressor(compressor);
        writer->Write();
        writer->Delete();

        compressor->Delete();
        dataSet->Delete();

        // Clear partition data
        X[part].clear();
        Y[part].clear();
        Z[part].clear();

        NDGLNO[part].clear();

        Ux[part].clear();
        Uy[part].clear();
        Uz[part].clear();

        data_arrays[part].clear();

        Pressure[part].clear();
      }
    }

    // Get out if we're finished
    if((args.find('m')==args.end())&&(npartitions>1))
      return 0;

    // Merge data if necessary
    if(args.find('m')!=args.end()){
      // Sanity checks
      if(npartitions<2){
        cerr<<"ERROR: Merge (-m) requested on a serial file. If you are operating on "
            <<"a serial file then drop the -m flag. Otherwise there may be a problem "
            <<"with your input.\n";
        exit(-1);
      }

      if(COLGA[0].size()==0){
        cerr<<"ERROR: This looks like an old fluidity output file - halo information "
            <<"is missing. Gerard refused to implement brute force merge on "
            <<"principle and went instead to find beer..\n";
        exit(-1);
      }

      // Global numbering offsets
      vector<size_t> offsets(npartitions);
      offsets[0] = 0;
      for(size_t i=1;i<npartitions;i++)
        offsets[i] = offsets[i-1] + n_private_nodes[i-1];

      // Initialise local-to-global lut
      deque< deque<int> > local_to_global(npartitions);
      for(size_t i=0; i<npartitions; i++){
        local_to_global[i].resize(X[i].size());
        fill(local_to_global[i].begin(), local_to_global[i].end(), -1);
        for(size_t j=0;j<n_private_nodes[i];j++){
          local_to_global[i][j] = offsets[i] + j;
        }
      }

      // Renumber COLGA
      for(size_t i=0;i<npartitions;i++){
        for(vector<int>::iterator it=COLGA[i].begin(); it!=COLGA[i].end(); it++){
          assert(*it>0);
          assert(*it<=(int)n_private_nodes[i]);
          *it = offsets[i] + *it - 1;
        }
      }

      // Correct local_to_global using halo information
      for(size_t i=0;i<npartitions;i++){
        for(size_t j=0;j<npartitions;j++){
          int nhalo=ATOREC[i][j+1]-ATOREC[i][j];
          assert(nhalo==(ATOSEN[j][i+1]-ATOSEN[j][i]));
          for(int k=0;k<nhalo;k++){
            int rindex = ATOREC[i][j] - 1 + k;
            int sindex = ATOSEN[j][i] - 1 + k;
            size_t nid = SCATER[i][rindex]-1;
            assert(nid>=n_private_nodes[i]);
            assert(nid<X[i].size());
            local_to_global[i][nid] = COLGA[j][sindex];
          }
        }
      }

      // Renumber and resize element lists
      size_t nloc = NDGLNO[0].size()/nelements[0];
      for(size_t i=0;i<npartitions;i++){

        // Establish who gets to write each element
        map<size_t, size_t> nowners;
        for(size_t j=0;j<npartitions;j++){
          int nhalo=ATOREC[i][j+1]-ATOREC[i][j];
          for(int k=0;k<nhalo;k++){
            int rindex = ATOREC[i][j] - 1 + k;
            size_t nid = SCATER[i][rindex]-1;
            assert(nowners.find(nid)==nowners.end());
            nowners[nid] = j;
          }
        }

        size_t elm_keep=0;
        for(size_t e=0;e<nelements[i];e++){
          int owner=npartitions+1;
          for(size_t j=0;j<nloc;j++){
            size_t nid = NDGLNO[i][e*nloc+j];
            if(nowners.find(nid)!=nowners.end()){
              owner = min(owner, (int)nowners[nid]);
            }else{
              owner = min(owner, (int)i);
            }
      if(local_to_global[i][NDGLNO[i][e*nloc+j]]<0){
        owner = -1;
      }
          }
    
          if(owner==i){
            for(size_t j=0;j<nloc;j++){
        NDGLNO[i][elm_keep*nloc+j] = local_to_global[i][NDGLNO[i][e*nloc+j]];
            }
            elm_keep++;
          }
        }
        NDGLNO[i].resize(elm_keep*nloc);
        nelements[i] = elm_keep;
      }

      // Merge arrays.
      size_t nphases = Ux[0].size();
      size_t ndata_arrays=data_arrays[0].size();

      // - first remove halo info
      for(size_t i=0; i<npartitions; i++){
        X[i].resize(n_private_nodes[i]);
        Y[i].resize(n_private_nodes[i]);
        Z[i].resize(n_private_nodes[i]);

        assert(nphases==Ux[i].size());
        assert(nphases==Uy[i].size());
        assert(nphases==Uz[i].size());

        for(size_t j=0; j<nphases; j++){
          Ux[i][j].resize(n_private_nodes[i]);
          Uy[i][j].resize(n_private_nodes[i]);
          Uz[i][j].resize(n_private_nodes[i]);
        }

        assert(ndata_arrays==data_arrays[i].size());
        for(size_t j=0; j<ndata_arrays; j++){
          data_arrays[i][j].resize(n_private_nodes[i]);
        }

        Pressure[i].resize(n_private_nodes[i]);
      }

      // - append all data to first partition, clearing data as we go.
      for(size_t i=1; i<npartitions; i++){
        X[0].insert(X[0].end(), X[i].begin(), X[i].end()); X[i].clear();
        Y[0].insert(Y[0].end(), Y[i].begin(), Y[i].end()); Y[i].clear();
        Z[0].insert(Z[0].end(), Z[i].begin(), Z[i].end()); Z[i].clear();

        assert(nphases==Ux[i].size());
        assert(nphases==Uy[0].size()); assert(nphases==Uy[i].size());
        assert(nphases==Uz[0].size()); assert(nphases==Uz[i].size());

        for(size_t j=0; j<nphases; j++){
          Ux[0][j].insert(Ux[0][j].end(), Ux[i][j].begin(), Ux[i][j].end()); Ux[i][j].clear();
          Uy[0][j].insert(Uy[0][j].end(), Uy[i][j].begin(), Uy[i][j].end()); Uy[i][j].clear();
          Uz[0][j].insert(Uz[0][j].end(), Uz[i][j].begin(), Uz[i][j].end()); Uz[i][j].clear();
        }
        Ux[i].clear();
        Uy[i].clear();
        Uz[i].clear();

        size_t ndata_arrays=data_arrays[0].size(); assert(ndata_arrays==data_arrays[i].size());
        for(size_t j=0; j<ndata_arrays; j++){
          data_arrays[0][j].insert(data_arrays[0][j].end(), data_arrays[i][j].begin(), data_arrays[i][j].end()); data_arrays[i][j].clear();
        }

        Pressure[0].insert(Pressure[0].end(), Pressure[i].begin(), Pressure[i].end()); Pressure[i].clear();
      }

      // Merge element lists
      for(size_t i=1;i<npartitions;i++){
        NDGLNO[0].insert(NDGLNO[0].end(), NDGLNO[i].begin(), NDGLNO[i].end());
        NDGLNO[i].clear();
      }

      // Clear these LUTS as they no longer mean anything when we've
      // merged. Periodic is simply not supported in this context.
      for(size_t i=0;i<npartitions;i++)
        periodic_renumbering[i].clear();
    }

    {
      if(verbose) cout<<"write out a serial vtu file\n";
      vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
      int nelm=nelements[0];
      for(size_t i=1;i<npartitions;i++)
        nelm+=nelements[i];
      Create_vtkUnstructuredGrid(dataSet, -1,
                                 Scalars, Vectors, Tensors, Tensors2,
                                 X.size(), nelm, vtk_cell_type,
                                 X[0], Y[0], Z[0],
                                 NDGLNO[0],
                                 Ux[0], Uy[0], Uz[0],
                                 data_arrays[0], Pressure[0],
                                 periodic_renumbering[0]);

      vtkXMLUnstructuredGridWriter *writer= vtkXMLUnstructuredGridWriter::New();
      vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
      char vtk_file_name[2048];
      sprintf(vtk_file_name,"%s_%d.vtu",project_name.c_str(),dump_id);
      writer->SetFileName(vtk_file_name);
      writer->SetInput(dataSet);
      writer->SetCompressor(compressor);
      writer->Write();
      writer->Delete();

      compressor->Delete();
      dataSet->Delete();
    }
  }

  return(0);
}
#else
int main(){
  std::cerr<<"No VTK support compiled with fluidity\n";
}
#endif
