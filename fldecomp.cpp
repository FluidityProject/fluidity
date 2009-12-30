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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <unistd.h>

#ifndef _AIX
#include <getopt.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "vtk.h"
#include "Halos_IO.h"
#include "Fldecomp_Wrappers.h"
#include "fmangle.h"
#include "partition.h"

using namespace std; 

using namespace Fluidity;

#define DDEBUG 1
#undef NDEBUG

extern "C" {
  void fldecomp_fc(const char *, const int *, const int *);
  void set_global_debug_level_fc(int *val);
}
 
void usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] -n nparts file\n"
      <<"\t-d\t\tPrint out partition diagnostics.\n"
      <<"\t-e <id>\t\tRecognise as extruded in the direction perpendicular "
      <<"to the surface with id.\n"
      <<"\t-f <file name>\t\tInput file (can alternatively specify as final "
      <<"argument)\n"
      <<"\t-h\t\tPrints out this message\n"
      <<"\t-i <file format>\t\tInput mesh type. Valid values are: triangle\n"
      <<"\t-k \t\tPartition a graph into k equal-size parts using the "
      <<"multilevel k-way partitioning algorithm (METIS PartGraphKway). This "
      <<"is the default if the number of partitions is greater than 8.\n"
      <<"\t-n <number of partitions>\t\tNumber of parts\n"
      <<"\t-r\t\tPartition a graph into k equal-size parts using multilevel recursive "
      <<"bisection (METIS PartGraphRecursive). This is the default if num partitions "
      <<"is less or equal to 8.\n"
      <<"\t-t\t\tRecognise as a Terreno output that is extruded in the Z direction\n"
      <<"\t-v\t\tVerbose mode\n";
  exit(-1);
}

vector<double> normal(const double *r0, const double *r1, const double *r2){
  double a = r1[0]-r0[0];
  double b = r1[1]-r0[1];
  double c = r1[2]-r0[2];

  double d = r2[0]-r0[0];
  double e = r2[1]-r0[1];
  double f = r2[2]-r0[2];

  vector<double> n(4);

  n[0] =  (b*f - e*c);
  n[1] = -(a*f - c*d);
  n[2] =  (a*e - b*d);

  double mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

  n[0] /= mag;
  n[1] /= mag;
  n[2] /= mag;
  n[3] = mag;

  return n;
}

void create_partitions(bool verbose,
                       const vector<int>& decomp, int nparts, int nnodes, int nelms, int nloc, const vector<int>& ENList,
                       vector<deque<int> >& nodes, vector<int>& npnodes, vector<deque<int> >& elements,
                       vector< vector<set<int> > >& halo1, vector< vector<set<int> > >& halo2){
  nodes.clear();    nodes.resize(nparts);
  npnodes.clear();  npnodes.resize(nparts);
  elements.clear(); elements.resize(nparts);

  halo1.clear();
  halo1.resize(nparts);
  for(int part=0; part<nparts; part++)
    halo1[part].resize(nparts);
  
  halo2.clear(); 
  halo2.resize(nparts);
  for(int part=0; part<nparts; part++)
    halo2[part].resize(nparts);
  
  for(int part=0; part<nparts; part++){
    if(verbose)
      cout<<"Making partition "<<part<<endl;
    
    // Find owned nodes 
    for(int nid=0; nid<nnodes; nid++){
      if(decomp[nid]==part)
        nodes[part].push_back(nid+1);
    }
    npnodes[part] = nodes[part].size();
    if(verbose)
      cout<<"Found "<<npnodes[part]<<" owned nodes\n";

    // Find elements with owned nodes and halo1
    set<int> halo_nodes;
    deque< pair<int, int> > sorted_elements;
    for(int eid=0;eid<nelms;eid++){
      int halo_count=0;
      pair<int, int> owned_elm(decomp[ENList[eid*nloc] - 1], eid);
      if(decomp[ENList[eid*nloc] - 1]!=part){
        halo_count++;
      }
      
      for(int j=1;j<nloc;j++){
        owned_elm.first = min(owned_elm.first, decomp[ENList[eid*nloc+j] - 1]);
        if(decomp[ENList[eid*nloc+j] - 1]!=part){
          halo_count++;
        }
      } 
      
      if(halo_count<nloc){
        sorted_elements.push_back(owned_elm);        
        if(halo_count>0){
          for(int j=0;j<nloc;j++){
            int nid = ENList[eid*nloc+j] - 1;
            if(decomp[nid]!=part){
              halo1[part][decomp[nid]].insert(nid+1);
              halo_nodes.insert(nid+1);
            }
          }
        }
      }
    }

    for(deque< pair<int, int> >::const_iterator it=sorted_elements.begin();it!=sorted_elements.end();++it)
      if(it->first==part)
        elements[part].push_back(it->second);
    for(deque< pair<int, int> >::const_iterator it=sorted_elements.begin();it!=sorted_elements.end();++it)
      if(it->first!=part)
        elements[part].push_back(it->second);

    if(verbose)
      cout<<"Found "<<halo_nodes.size()<<" halo1 nodes\n";
    
    // Find halo2 elements and nodes
    set<int> halo2_elements,  more_halo_nodes;
    for(int eid=0; eid<nelms; eid++){
      int owned_node_count=0;
      bool touches_halo1=false;
      for(int j=0;j<nloc;j++){
        int fnid = ENList[eid*nloc+j];
        if(halo_nodes.count(fnid)){
          touches_halo1 = true;
        }
        if(decomp[fnid-1]==part)
          owned_node_count++;
      }

      if(touches_halo1&&(owned_node_count==0)){
        halo2_elements.insert(eid);
        for(int j=0;j<nloc;j++){
          int fnid = ENList[eid*nloc+j];
          if(!halo_nodes.count(fnid)){
            halo2[part][decomp[fnid-1]].insert(fnid);
            more_halo_nodes.insert(fnid);
          }
        }
      }
    }

    if(verbose)
      cout<<"Found "<<halo2_elements.size()<<" halo2 elements\n";

    for(int i=0;i<nparts;i++)
      for(set<int>::const_iterator it=halo1[part][i].begin(); it!=halo1[part][i].end(); ++it)
        halo2[part][i].insert(*it);

    for(set<int>::const_iterator it=halo_nodes.begin(); it!=halo_nodes.end(); ++it)
      nodes[part].push_back(*it);
    
    for(set<int>::const_iterator it=more_halo_nodes.begin(); it!=more_halo_nodes.end(); ++it)
      nodes[part].push_back(*it);
    
    for(set<int>::const_iterator it=halo2_elements.begin(); it!=halo2_elements.end(); ++it)
      elements[part].push_back(*it);

    if(verbose)
      cout<<"Partition: "<<part<<", Private nodes: "<<npnodes[part]<<", Total nodes: "<<nodes[part].size()<<"\n";
  }
  
  return;
}

int main(int argc, char **argv){
 #ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI::Init(argc,argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
 
  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> flArgs;
  while((c = getopt(argc, argv, "de:f:hi:kn:rtv")) != -1){
    if (c != '?'){
      if (optarg == NULL){
        flArgs[c] = "true";
      }else{
        flArgs[c] = optarg;
      }
    }else{
      if (isprint(optopt)){
        cerr << "Unknown option " << optopt << endl;
      }else{
        cerr << "Unknown option " << hex << optopt << endl;
      }
      usage(argv[0]);
      exit(-1);
    }
  }

  // Help?
  if(flArgs.count('h')||(flArgs.count('n')==0)){
    usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  bool verbose=false;
  int val=3;
  if(flArgs.count('v')==0)
    val = 0;
  else
    verbose = true;
  set_global_debug_level_fc(&val);
  
  if(flArgs.find('f')==flArgs.end()){
    if(argc>optind+1){
      flArgs['f'] = argv[optind+1];
    }else if(argc==optind+1){
      flArgs['f'] = argv[optind];
    }
  }
  
  string filename = flArgs['f'];

  string file_format;
  if(flArgs.find('i')==flArgs.end()){
     cerr<<"Warning: No mesh type specified - please supply a mesh type"<<endl;
     if(filename.substr(filename.size()-4, 4)==string(".dat")){
       file_format="gid_dat";
     }else{
       file_format="triangle";
     }
  }else{
    file_format=flArgs['i'];
  }
  if(verbose)
    cout<<"Input mesh type: "<<file_format<<endl;
  
  vector<int> surface_nids;
  int nparts = atoi(flArgs['n'].c_str());
  if(nparts<2){
    cerr<<"ERROR: number of partitions requested is less than 2. Please check your usage and try again.\n";
    usage(argv[0]);
    exit(-1);
  }

  if(file_format=="triangle"){
    const int halo1_level = 1, halo2_level = 2;
    
    // Read in the mesh
    if(verbose)
      cout<<"Reading in triangle mesh with base name "<<filename<<"\n";
    vector<double> x;
    int dim;
    vector<int> ENList, regionIds;
    int nloc;
    deque<vector<int> > SENList;
    vector<int> boundaryIds;
    int snloc;
    if(ReadMesh(filename, file_format, x, dim,
                ENList, regionIds, nloc, SENList, boundaryIds, snloc)){
      cerr<<"ERROR: failed to read mesh file with base name "<<filename<<endl;
      exit(-1);
    }
    
    int nnodes=x.size()/dim;
    vector<int> surface_nids;
    
    if(flArgs.find('e')!=flArgs.end()){
      // This triangle file has been extruded perpendicular for
      // surface with the id flArgs.find('e'). The domain
      // decomposition goes along this surface.
      
      // Initialise surface_nids
      surface_nids.resize(nnodes);
      for(int i=0;i<nnodes;i++)
        surface_nids[i] = -1;
      
      int eplane = atoi(flArgs['e'].c_str());
      map<int, set<int> > sNSList;
      map<int, vector<double> > eplane_elements;
      int nselements = SENList.size();
      for(int i=0;i<nselements;i++){
        if(boundaryIds[i]==eplane){
          eplane_elements[i] = normal(&(x[3*(SENList[i][0]-1)]),
                                      &(x[3*(SENList[i][1]-1)]),
                                      &(x[3*(SENList[i][2]-1)]));
          for(int j=0;j<snloc;j++){
            sNSList[SENList[i][j]-1].insert(i);
            surface_nids[SENList[i][j]-1] = SENList[i][j]-1;
          }
        }
      }

      // Approximate the normal at each vertex
      map<int, vector<double> > normals;
      for(map<int, set<int> >::const_iterator it=sNSList.begin();it!=sNSList.end();++it){
        double patch_area = 0.0;
        vector<double> weighed_normal(3, 0.0);
        for(set<int>::const_iterator jt=it->second.begin();jt!=it->second.end();++jt){
          patch_area += eplane_elements[*jt][3];
          for(int i=0;i<3;i++)
            weighed_normal[i] += eplane_elements[*jt][i]*eplane_elements[*jt][3];
        }
        
        // renormalise
        double mag = sqrt(weighed_normal[0]*weighed_normal[0]+
                          weighed_normal[1]*weighed_normal[1]+
                          weighed_normal[2]*weighed_normal[2]);
        for(int i=0;i<3;i++)
          weighed_normal[i] /= mag;
        
        normals[it->first] = weighed_normal;
      }
      
      // Get the useful node-node adjancy list
      vector< set<int> > NNList(nnodes);
      for(int i=0;i<nnodes;i++){
        for(int j=0;j<nloc;j++){
          if(surface_nids[ENList[i*nloc+j]-1]==-1)
            NNList[i].insert(ENList[i*nloc+j]-1);
        }
      }
      
      // Finally, discover the extrusion.
      for(map<int, set<int> >::const_iterator it=sNSList.begin();it!=sNSList.end();++it){
        int enode = it->first;
        set<int>::const_iterator next_node=NNList[enode].begin();
        for(;next_node!=NNList[enode].end();++next_node){
          double dot = 0;
          for(int i=0;i<3;i++)
            dot += normals[enode][i]*(x[enode*3+i]-x[(*next_node)*3+i]);
          
          if(fabs(fabs(dot)-1)<0.0e-5){
            break;
          }
        }
        if(next_node==NNList[enode].end()){
          break;
        }else{
          surface_nids[*next_node] = surface_nids[enode];
          enode = *next_node;
        }
      }

      // Check all is well
      int err=0;
      for(int i=0;i<nnodes;i++){
        if(surface_nids[i]<0){
          err=-1;
          cerr<<"ERROR: failed to cover one of the extruded nodes.\n";
        }
      }
      if(err)
        exit(-1);
    }else if(flArgs.find('t')!=flArgs.end()){
      // This triangle file came from Terreno and should have
      // attribure data indicating the surface node lieing above the
      // node in the extruded mesh.
      string filename_node = filename+".node";
      fstream afile;
      afile.open(filename_node.c_str(), ios::in);
      if(!afile.is_open()){
        cerr<<"ERROR: Shoreline data file, "<<filename_node
            <<", cannot be opened. Does it exist? Have you read permission?\n";
        exit(-1);
      }
      afile.close();
      
      FILE *file=fopen(filename_node.c_str(), "r");
      
      int header[4];
      for(int i=0;i<4;i++)
        fscanf(file, "%d", header+i);
      if(header[2]==1){
        assert(nnodes==header[0]);
        surface_nids.resize(nnodes);
        double dummy;
        for(int i=0;i<nnodes;i++){
          for(int j=0;j<=header[1];j++)
            fscanf(file, "%lf", &dummy);
          fscanf(file, "%d", &(surface_nids[i]));
        }
      }
      fclose(file);
    }
    
    vector<int> decomp;
    int partition_method;
    if(nparts<=8)
      partition_method = 0; // METIS PartGraphRecursive
    else
      partition_method = 1; // METIS PartGraphKway

    if(flArgs.find('r')!=flArgs.end()){
      partition_method = 0; // METIS PartGraphRecursive
      
      if(flArgs.find('k')!=flArgs.end()){
        cerr<<"WARNING: should not specify both -k and -r. Choosing -r.\n";
      }
    }
    if(flArgs.find('k')!=flArgs.end()){
      partition_method = 1; // METIS PartGraphKway
    }

    int edgecut=0;
    if(surface_nids.size()){
      // Partition the mesh
      if(verbose)
        cout<<"Partitioning the extruded Terreno mesh\n";

      // Partition the mesh. Generates a map "decomp" from node number
      // (numbered from zero) to partition number (numbered from
      // zero).
      edgecut = partition(ENList, surface_nids, nloc, nnodes, nparts, partition_method, decomp);
    }else{
      // Partition the mesh
      if(verbose)
        cout<<"Partitioning the mesh\n";
      
      // Partition the mesh. Generates a map "decomp" from node number
      // (numbered from zero) to partition number (numbered from
      // zero).
      edgecut = partition(ENList, dim, nloc, nnodes, nparts, partition_method, decomp);
    }

    if(flArgs.find('d')!=flArgs.end()){
      cout<<"Edge-cut: "<<edgecut<<endl;
    }

    // Process the partitioning
    if(verbose)
      cout<<"Processing the mesh partitions\n";
    vector<deque<int> > nodes(nparts);
    vector<int> npnodes(nparts);
    vector<deque<int> > elements(nparts);
    vector< vector<set<int> > > halo1(nparts), halo2(nparts);
    // Construct:
    //   nodes    - nodes in each partition (private and halo), numbered from
    //              one
    //   npnodes  - number of nodes private to each partition
    //   elements - elements ids in each partition (private and halo), numbered
    //              from zero
    //   halo1    - receive nodes (numbered from one) in the first halo for each
    //              partition
    //   halo2    - receive nodes (numbered from one) in the second halo for
    //              each partition
    // where in each case partitions are numbered from zero.
    create_partitions(verbose,
                      decomp, nparts, nnodes, ENList.size()/nloc, nloc, ENList,
                      nodes, npnodes, elements, halo1, halo2);
                 
    // Extract and write out data for each partition
    for(int i=0;i<nparts;i++){
      if(verbose)
        cout<<"Extracting mesh data for partition "<<i<<"\n";
      
      // Map from global node numbering (numbered from one) to partition node 
      // numbering (numbered from one)
      map<int, int> renumber;
      for(size_t j=0;j<nodes[i].size();j++){
        assert(renumber.find(nodes[i][j])==renumber.end());
        renumber[nodes[i][j]]=j+1;
      }
      
      // Coordinate data
      vector<double> partX(nodes[i].size()*dim);
      for(size_t j=0;j<nodes[i].size();j++){
        for(int k=0;k<dim;k++){
          partX[j * dim + k] = x[(nodes[i][j] - 1) * dim + k];
        }
      }
      
      // Volume element data
      vector<int> partENList;
      vector<int> partRegionIds;
      // Map from partition node numbers (numbered from one) to partition
      // element numbers (numbered from zero)
      multimap<int, int> partNodesToEid;
      for(deque<int>::const_iterator iter=elements[i].begin();iter!=elements[i].end();iter++){
        for(int j=0;j<nloc;j++){
          assert(renumber.find(ENList[*iter*nloc+j])!=renumber.end());
          partENList.push_back(renumber[ENList[*iter*nloc+j]]);
          partNodesToEid.insert(pair<int, int>(*(partENList.rbegin()),(partENList.size()-1)/nloc));
        }
        partRegionIds.push_back(regionIds[*iter]);
      }
      
      // Surface element data
      deque<vector<int> > partSENList;
      vector<int> partBoundaryIds;
#ifndef NDEBUG
      set< vector<int> > surf_elms;
#endif
      for(size_t j=0;j<SENList.size();j++){
        // In order for a global surface element to be a partition surface
        // element, all of its nodes must be attached to at least one partition
        // volume element
        if(SENList[j].size()==0 or renumber.find(SENList[j][0])==renumber.end() or partNodesToEid.count(renumber[SENList[j][0]])==0){
          continue;
        }
        
        bool SEOwned=false;
        for(multimap<int, int>::const_iterator iter=partNodesToEid.find(renumber[SENList[j][0]]);iter->first==renumber[SENList[j][0]];iter++){
          SEOwned=true;
          set<int> VENodes;
          for(int k=iter->second*nloc;k<iter->second*nloc+nloc;k++){
            VENodes.insert(partENList[k]);
          }
          for(size_t k=1;k<SENList[j].size();k++){
            if(renumber.find(SENList[j][k])==renumber.end() or VENodes.count(renumber[SENList[j][k]])==0){
              SEOwned=false;
              break;
            }
          }
          if(SEOwned){
            break;
          }
        }

        if(SEOwned){
          vector<int> elm(SENList[j].size());
          
          for(size_t k=0;k<SENList[j].size();k++){
            assert(renumber.find(SENList[j][k])!=renumber.end());
            elm[k]=renumber[SENList[j][k]];
          }
#ifndef NDEBUG
          assert(surf_elms.find(elm)==surf_elms.end());
          surf_elms.insert(elm);
#endif
          partSENList.push_back(elm);
          partBoundaryIds.push_back(boundaryIds[j]);
        }
      }
      
      // Write out the partition mesh
      ostringstream buffer;
      buffer<<filename<<"_"<<i;
      if(verbose)
        cout<<"Writing out triangle mesh for partition "<<i<<" to files with base name "<<buffer.str()<<"\n";
      if(WriteMesh(buffer.str(), file_format, partX, dim,
                   partENList, partRegionIds, nloc, partSENList, partBoundaryIds, snloc)){
        cerr<<"ERROR: failed to write mesh file with base name "<<buffer.str()<<endl;
        exit(-1);
      }
      buffer.str("");

      // Extract halo data
      if(verbose)
        cout<<"Extracting halo data for partition "<<i<<"\n";
      map<int, vector< vector<int> > > send, recv;
      map<int, int> npnodes_handle;

      recv[halo1_level].resize(nparts);
      send[halo1_level].resize(nparts);

      recv[halo2_level].resize(nparts);
      send[halo2_level].resize(nparts);

      for(int j=0;j<nparts;j++){
        for(set<int>::const_iterator it=halo1[i][j].begin();it!=halo1[i][j].end();++it){
          assert(renumber.find(*it)!=renumber.end());
          assert(renumber[*it]>npnodes[i]);
          recv[halo1_level][j].push_back(renumber[*it]);
        }
        for(set<int>::const_iterator it=halo1[j][i].begin();it!=halo1[j][i].end();++it){
          assert(renumber.find(*it)!=renumber.end());
          assert(renumber[*it]<=npnodes[i]);
          send[halo1_level][j].push_back(renumber[*it]);
        }

        for(set<int>::const_iterator it=halo2[i][j].begin();it!=halo2[i][j].end();++it){
          assert(renumber.find(*it)!=renumber.end());
          assert(renumber[*it]>npnodes[i]);
          recv[halo2_level][j].push_back(renumber[*it]);
        }
        for(set<int>::const_iterator it=halo2[j][i].begin();it!=halo2[j][i].end();++it){
          assert(renumber.find(*it)!=renumber.end());
          assert(renumber[*it]<=npnodes[i]);
          send[halo2_level][j].push_back(renumber[*it]);
        }
      }

      npnodes_handle[halo1_level]=npnodes[i];
      npnodes_handle[halo2_level]=npnodes[i];

      buffer<<filename<<"_"<<i<<".halo";
      if(verbose)
        cout<<"Writing out halos for partition "<<i<<" to file "<<buffer.str()<<"\n";

      if(WriteHalos(buffer.str(), i, nparts, npnodes_handle, send, recv)){
        cerr<<"ERROR: failed to write halos to file "<<buffer.str()<<endl;
        exit(-1);
      }
      buffer.str("");
    }
  }else{
    cerr<<"ERROR: unrecognised file type\n";
  }
  
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return(0);
}
