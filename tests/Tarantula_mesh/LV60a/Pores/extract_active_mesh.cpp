#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <limits>

#include <cassert>
#include <cstdlib>

#include <getopt.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

void usage(char *cmd){
  std::cerr<<"\nUsage: "<<cmd<<" [options ...] [Tarantula mesh file]\n"
           <<"\nOptions:\n"
           <<" -h, --help\n\tHelp! Prints this message.\n"
           <<" -v, --verbose\n\tVerbose output.\n"
           <<" -b ID0,ID1, --boundary ID0, ID1\n\tBoundary ID's corresponding to the inflow and outflow.\n";
  return;
}

int parse_arguments(int argc, char **argv,
                    std::string &filename, bool &verbose, int *IDs){

  // Set defaults
  verbose = false;
  IDs[0]=-1;
  IDs[0]=-1;

  if(argc==1){
    usage(argv[0]);
    exit(0);
  }

  struct option longOptions[] = {
    {"help", 0, 0, 'h'},
    {"verbose", 0, 0, 'v'},
    {"boundary", optional_argument, 0, 'b'},
    
    {0, 0, 0, 0}
  };

  int optionIndex = 0;
  int verbosity = 0;
  int c;
  const char *shortopts = "hvb:";

  std::string boundary_str;
  size_t pos;
  std::stringstream id0, id1;

  // Set opterr to nonzero to make getopt print error messages
  opterr=1;
  while (true){
    c = getopt_long(argc, argv, shortopts, longOptions, &optionIndex);
    
    if (c == -1) break;
    
    switch (c){
    case 'h':
      usage(argv[0]);
      break;
    case 'v':
      verbose = true;
      break;
    case 'b':
      boundary_str = std::string(optarg);
      pos = boundary_str.find(",", 0);
      if(pos == std::string::npos){
        std::cerr<<"ERROR: boundary not specified correctly.\n";
        usage(argv[0]);
        exit(0);
      }
      id0<<boundary_str.substr(0, pos);
      if(!(id0>>IDs[0])) // Give the value to Result using the characters in the string
        IDs[0] = -1;
      id1<<boundary_str.substr(pos+1);
      if(!(id1>>IDs[1])) // Give the value to Result using the characters in the string
        IDs[1] = -1;

      // parse optarg
      break;    
    case '?':
      // missing argument only returns ':' if the option string starts with ':'
      // but this seems to stop the printing of error messages by getopt?
      std::cerr<<"ERROR: unknown option or missing argument\n";
      usage(argv[0]);
      exit(-1);
    case ':':
      std::cerr<<"ERROR: missing argument\n";
      usage(argv[0]);
      exit(-1);
    default:
      // unexpected:
      std::cerr<<"ERROR: getopt returned unrecognized character code\n";
      exit(-1);
    }
  }

  filename = std::string(argv[argc-1]);
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.good()){
    std::cerr<<"ERROR: Cannot read file: "<<filename<<std::endl;
    usage(argv[0]);
    exit(0);
  }

}

int read_tarantula_mesh_file(std::string filename,
                             std::vector<double> &xyz,
                             std::vector<int> &tets, 
                             std::vector<int> &facets,
                             std::vector<int> &facet_ids){
  std::ifstream infile;
  infile.open(filename.c_str());
  
  // Read header
  std::string throwaway;
  std::getline(infile, throwaway); // line 1
  std::getline(infile, throwaway); // line 2
  int NNodes;
  infile>>NNodes;
  
  // Read vertices
  xyz.resize(NNodes*3);
  for(int i=0;i<NNodes;i++){
    infile>>xyz[i*3];
    infile>>xyz[i*3+1]; 
    infile>>xyz[i*3+2];
  }

  // throwaway trash
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway); 
  std::getline(infile, throwaway);

  // Read elements
  int NTetra, nloc;
  infile>>NTetra;
  tets.resize(NTetra*4);
  for(int i=0;i<NTetra;i++){
    infile>>nloc;
    assert(nloc==4);
    infile>>tets[i*4];
    infile>>tets[i*4+1];
    infile>>tets[i*4+2];
    infile>>tets[i*4+3];
  }

  // Read facets
  while(infile.good()){
    // Stream through file until we find facet data.
    std::getline(infile, throwaway); 
    if(throwaway.substr(0, 4)!="Face")
      continue;

    // Get facet ID
    std::stringstream id_ss(throwaway.substr(4));
    int id;
    id_ss>>id;

    std::getline(infile, throwaway); 
    int nfacets;
    infile>>nfacets;
    nfacets/=2;

    for(int i=0;i<nfacets;i++){
      int eid, index;
      infile>>eid;
      infile>>index;
      assert(index<4);

      facet_ids.push_back(id);
      if(index==0){
        facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+2]);
      }else if(index==1){
        facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+3]);
      }else if(index==2){
        facets.push_back(tets[eid*4+2]); facets.push_back(tets[eid*4+1]); facets.push_back(tets[eid*4+3]);
      }else if(index==3){
        facets.push_back(tets[eid*4+0]); facets.push_back(tets[eid*4+2]); facets.push_back(tets[eid*4+3]);
      }
    }
  }

  infile.close();
}

int write_vtk_file(std::string filename,
                   std::vector<double> &xyz,
                   std::vector<int> &tets, 
                   std::vector<int> &facets,
                   std::vector<int> &facet_ids){
  
  // Write out points
  int NNodes = xyz.size()/3;
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NNodes);
  for(int i=0;i<NNodes;i++)
    pts->SetPoint(i, &(xyz[i*3]));
  
  // Initalise the vtk mesh
  vtkUnstructuredGrid *ug_tets = vtkUnstructuredGrid::New();
  ug_tets->SetPoints(pts);

  int NTetra = tets.size()/4;
  for(int i=0;i<NTetra;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<4;j++)
      idlist->InsertNextId(tets[i*4+j]);
    ug_tets->InsertNextCell(10, idlist);
    idlist->Delete();
  }
  
  vtkXMLUnstructuredGridWriter *tet_writer = vtkXMLUnstructuredGridWriter::New();
  tet_writer->SetFileName(std::string(filename+".vtu").c_str());
  tet_writer->SetInput(ug_tets);
  tet_writer->Write();
  
  ug_tets->Delete();
  tet_writer->Delete();

  // Write out facets
  vtkUnstructuredGrid *ug_facets = vtkUnstructuredGrid::New();
  ug_facets->SetPoints(pts);
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    vtkIdList *idlist = vtkIdList::New();
    for(int j=0;j<3;j++){
      idlist->InsertNextId(facets[i*3+j]);
    }
    ug_facets->InsertNextCell(5, idlist);
    idlist->Delete();
  }

  vtkIntArray *vtk_facet_ids = vtkIntArray::New();
  vtk_facet_ids->SetNumberOfTuples(NFacets);
  vtk_facet_ids->SetNumberOfComponents(1);
  vtk_facet_ids->SetName("Facet IDs");
  for(int i=0;i<NFacets;i++){
    vtk_facet_ids->SetValue(i, facet_ids[i]);
  }
  ug_facets->GetCellData()->AddArray(vtk_facet_ids);
  vtk_facet_ids->Delete();
  
  vtkXMLUnstructuredGridWriter *tri_writer = vtkXMLUnstructuredGridWriter::New();
  tri_writer->SetFileName(std::string(filename+"_facets.vtu").c_str());
  tri_writer->SetInput(ug_facets);
  tri_writer->Write();

  pts->Delete();
  ug_facets->Delete();
  tri_writer->Delete();

  return 0;
}

int trim_channels(const int *id,
                  std::vector<double> &xyz,
                  std::vector<int> &tets, 
                  std::vector<int> &facets,
                  std::vector<int> &facet_ids){

  // Create node-element adjancy list.
  std::vector< std::set<int> > NEList(xyz.size()/3);
  int NTetra = tets.size()/4;
  for(int i=0;i<NTetra;i++){
    for(int j=0;j<4;j++){
      NEList[tets[i*4+j]].insert(i);
    }
  }
  
  // Create element-element adjancy list.
  std::vector<int> EEList(NTetra*4, -1);
  for(int i=0;i<NTetra;i++){
    for(int j=0;j<4;j++){
      
      std::set<int> edge_neighbours;
      set_intersection(NEList[tets[i*4+(j+1)%4]].begin(), NEList[tets[i*4+(j+1)%4]].end(),
                       NEList[tets[i*4+(j+2)%4]].begin(), NEList[tets[i*4+(j+2)%4]].end(),
                       inserter(edge_neighbours, edge_neighbours.begin()));
      
      std::set<int> neighbours;
      set_intersection(NEList[tets[i*4+(j+3)%4]].begin(), NEList[tets[i*4+(j+3)%4]].end(),
                       edge_neighbours.begin(), edge_neighbours.end(),
                       inserter(neighbours, neighbours.begin()));

      if(neighbours.size()==2){
        if(*neighbours.begin()==i)
          EEList[i*4+j] = *neighbours.rbegin();
        else
          EEList[i*4+j] = *neighbours.begin();
      }
    }
  }

  // Create full facet ID list. Also, create the initial fronts for
  // the active region detection.
  std::map< std::set<int>, int> facet_id_lut;
  int NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    std::set<int> facet;
    for(int j=0;j<3;j++){
      facet.insert(facets[i*3+j]);
    }
    assert(facet_id_lut.find(facet)==facet_id_lut.end());
    facet_id_lut[facet] = facet_ids[i];
  }

  std::set<int> front0, front1;
  std::vector<int> full_facet_id_list(NTetra*4, -1);
  for(int i=0;i<NTetra;i++){
    for(int j=0;j<4;j++){
      if(EEList[i*4+j]==-1){
        std::set<int> facet;
        for(int k=1;k<4;k++)
          facet.insert(tets[i*4+(j+k)%4]);
        
        std::map< std::set<int>, int>::iterator facet_id_pair = facet_id_lut.find(facet);
        if(facet_id_pair==facet_id_lut.end()){
          full_facet_id_list[i*4+j] = 0;
          for(int k=1;k<4;k++)
            facets.push_back(tets[i*4+(j+k)%4]);
          facet_ids.push_back(0);
        }else{
          if(facet_id_pair->second==id[0])
            front0.insert(i);
          else if(facet_id_pair->second==id[1])
            front1.insert(i);
          full_facet_id_list[i*4+j] = facet_id_pair->second;
        }
      }
    }
  }
  
  // Advance front0
  std::vector<int> label(NTetra, 0);
  while(!front0.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front0.begin();
    front0.erase(front0.begin());
    if(label[seed]==1)
      continue;
    label[seed] = 1;
    
    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]!=1){
        front0.insert(eid);
      }
    }
  }

  // Advance backsweep using front1.
  while(!front1.empty()){
    // Get the next unprocessed element in the set.
    int seed = *front1.begin();
    front1.erase(front1.begin());
    if(label[seed]!=1) // ie was either never of interest or has been processed in the backsweep.
      continue;
    label[seed] = 2;
    
    for(int i=0;i<4;i++){
      int eid = EEList[seed*4+i];
      if(eid!=-1 && label[eid]==1){
        front1.insert(eid);
      }
    }
  }

  // Find active vertex set and create renumbering.
  std::map<int, int> renumbering;
  for(int i=0;i<NTetra;i++){
    if(label[i]==2){
      for(int j=0;j<4;j++)
        renumbering.insert(std::pair<int, int>(tets[i*4+j], -1));
    }
  }

  // Create new compressed mesh.
  std::vector<double> xyz_new;
  std::vector<int> tets_new;
  std::vector<int> facets_new;
  std::vector<int> facet_ids_new;
  int cnt=0;
  for(std::map<int, int>::iterator it=renumbering.begin();it!=renumbering.end();++it){
    it->second = cnt++;
    
    xyz_new.push_back(xyz[(it->first)*3]);
    xyz_new.push_back(xyz[(it->first)*3+1]);
    xyz_new.push_back(xyz[(it->first)*3+2]);
  }
  for(int i=0;i<NTetra;i++){
    if(label[i]==2){
      for(int j=0;j<4;j++){
        tets_new.push_back(renumbering[tets[i*4+j]]);
      }
    }
  }
  NFacets = facet_ids.size();
  for(int i=0;i<NFacets;i++){
    int facet[3];
    bool redundant=false;
    for(int j=0;j<3;j++){
      std::map<int, int>::iterator it=renumbering.find(facets[i*3+j]);
      if(it!=renumbering.end()){
        facet[j] = it->second;
      }else{
        redundant = true;
        break;
      }
    }
    if(redundant)
      continue;
    for(int j=0;j<3;j++)
      facets_new.push_back(facet[j]);
    facet_ids_new.push_back(facet_ids[i]);
  }
  xyz.swap(xyz_new);
  tets.swap(tets_new);
  facets.swap(facets_new);
  facet_ids.swap(facet_ids_new);
}

int write_triangle_file(std::string basename,
                        std::vector<double> &xyz,
                        std::vector<int> &tets, 
                        std::vector<int> &facets,
                        std::vector<int> &facet_ids){
  std::string filename_node = basename+".node";
  std::string filename_face = basename+".face";
  std::string filename_ele = basename+".ele";
  
  int NNodes = xyz.size()/3;
  int NTetra = tets.size()/4;
  int NFacets = facet_ids.size();
  
  ofstream nodefile;
  nodefile.open(std::string(basename+".node").c_str());
  nodefile<<NNodes<<" "<<3<<" "<<0<<" "<<0<<std::endl;
  nodefile<<std::setprecision(std::numeric_limits<double>::digits10+1);

  for(int i=0;i<NNodes;i++){
    nodefile<<i+1<<" "<<xyz[i*3]<<" "<<xyz[i*3+1]<<" "<<xyz[i*3+2]<<std::endl;
  }

  ofstream elefile;
  elefile.open(std::string(basename+".ele").c_str());
  elefile<<NTetra<<" "<<4<<" "<<1<<std::endl;

  for(int i=0;i<NTetra;i++){
    elefile<<i+1<<" "<<tets[i*4]+1<<" "<<tets[i*4+1]+1<<" "<<tets[i*4+2]+1<<" "<<tets[i*4+3]+1<<" 1"<<std::endl;
  }

  ofstream facefile;
  facefile.open(std::string(basename+".face").c_str());
  facefile<<NFacets<<" "<<1<<std::endl;
  for(int i=0;i<NFacets;i++){
    facefile<<i+1<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<" "<<facet_ids[i]<<std::endl;
  }


  
  
  return 0;
}


int main(int argc, char **argv){
  std::string filename;
  bool verbose;
  int IDs[2];
  parse_arguments(argc, argv, filename, verbose, IDs);

  std::vector<double> xyz;
  std::vector<int> tets, facets, facet_ids;
  std::string basename = filename.substr(0, filename.size()-4);
  std::cout<<"basename = "<<basename<<std::endl;

  if(verbose) std::cout<<"INFO: Reading "<<filename<<std::endl;
  read_tarantula_mesh_file(filename, xyz, tets, facets, facet_ids);
  if(verbose) std::cout<<"INFO: Finished reading "<<filename<<std::endl;

  write_vtk_file(basename+"_original", xyz, tets, facets, facet_ids);

  if(verbose) std::cout<<"INFO: Trimming inactive regions."<<filename<<std::endl;
  trim_channels(IDs, xyz, tets, facets, facet_ids);
  if(verbose) std::cout<<"INFO: Finished trimming."<<filename<<std::endl;

  write_vtk_file(basename, xyz, tets, facets, facet_ids);
  write_triangle_file(basename, xyz, tets, facets, facet_ids);

  return 0;
}
