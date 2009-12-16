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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <set>

#include "fmangle.h"

using namespace std;

extern "C" { // This is the glue between METIS and fluidity  
  
  // Declarations needed from METIS
  typedef int idxtype;
  void METIS_PartGraphKway(int *,idxtype *,idxtype *,idxtype *,idxtype *,int *,int *,int *,
                           int *,int *,idxtype *);
  void METIS_PartGraphRecursive(int *,idxtype *,idxtype *,idxtype *,idxtype *,int *,int *,
                                int *,int *,int *,idxtype *);

#ifndef HAVE_MPI
        void METIS_PartGraphKway(int *a, idxtype *b, idxtype *c, idxtype *d, idxtype *e,
                                 int *f, int *g, int *h, int *i, int *j, idxtype *k){}
        void METIS_PartGraphRecursive(int *a, idxtype *b, idxtype *c, idxtype *d,
                                      idxtype *e, int *f, int *g, int *h, int *i, int *j, idxtype *k){}
#endif
}

namespace Fluidity{

  int FormGraph(const vector<int>& ENList, const int& dim, const int& nloc, const int& nnodes,
                vector<set<int> >& graph){
    int num_elems = ENList.size()/nloc;
  
    graph.clear();  graph.resize(nnodes);
    
    switch (dim){
    case 1:
      switch(nloc){
      case 2:
        // Lines
        for(int i=0;i<num_elems;i++){
          graph[ENList[i*nloc]-1].insert(ENList[i*nloc+1]);
          graph[ENList[i*nloc+1]-1].insert(ENList[i*nloc]);
        }
        break;
      default:
        cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
        return -1;
      }
      break;
    case 2:
      switch(nloc){
      case 3:
        // Linear triangles
        for(int i=0;i<num_elems;i++){
          for(int j=0;j<nloc;j++){
            for(int k=0;k<nloc;k++){
              if(j!=k)
                graph[ENList[i*nloc+j]-1].insert(ENList[i*nloc+k]);
            }
          }
        }
        break;
      default:
        cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
        return -1;
      }
      break;
    case 3:
      switch (nloc){
      case 4:
        // Linear tets
        for(int i=0;i<num_elems;i++){
          for(int j=0;j<nloc;j++){
            for(int k=0;k<nloc;k++){
              if(j!=k)
                graph[ENList[i*nloc+j]-1].insert(ENList[i*nloc+k]);
            }
          }
        }
        break;
      case 8:
        // Linear hexes
        for(int i=0;i<num_elems;i++){
          graph[ENList[i*nloc  ]-1].insert(ENList[i*nloc+1]);
          graph[ENList[i*nloc  ]-1].insert(ENList[i*nloc+2]);
          graph[ENList[i*nloc  ]-1].insert(ENList[i*nloc+4]);
          
          graph[ENList[i*nloc+1]-1].insert(ENList[i*nloc  ]);
          graph[ENList[i*nloc+1]-1].insert(ENList[i*nloc+3]);
          graph[ENList[i*nloc+1]-1].insert(ENList[i*nloc+5]);
          
          graph[ENList[i*nloc+2]-1].insert(ENList[i*nloc  ]);
          graph[ENList[i*nloc+2]-1].insert(ENList[i*nloc+3]);
          graph[ENList[i*nloc+2]-1].insert(ENList[i*nloc+6]);
          
          graph[ENList[i*nloc+3]-1].insert(ENList[i*nloc+1]);
          graph[ENList[i*nloc+3]-1].insert(ENList[i*nloc+2]);
          graph[ENList[i*nloc+3]-1].insert(ENList[i*nloc+7]);
          
          graph[ENList[i*nloc+4]-1].insert(ENList[i*nloc  ]);
          graph[ENList[i*nloc+4]-1].insert(ENList[i*nloc+5]);
          graph[ENList[i*nloc+4]-1].insert(ENList[i*nloc+6]);
          
          graph[ENList[i*nloc+5]-1].insert(ENList[i*nloc+1]);
          graph[ENList[i*nloc+5]-1].insert(ENList[i*nloc+4]);
          graph[ENList[i*nloc+5]-1].insert(ENList[i*nloc+7]);
          
          graph[ENList[i*nloc+6]-1].insert(ENList[i*nloc+2]);
          graph[ENList[i*nloc+6]-1].insert(ENList[i*nloc+4]);
          graph[ENList[i*nloc+6]-1].insert(ENList[i*nloc+7]);
          
          graph[ENList[i*nloc+7]-1].insert(ENList[i*nloc+3]);
          graph[ENList[i*nloc+7]-1].insert(ENList[i*nloc+5]);
          graph[ENList[i*nloc+7]-1].insert(ENList[i*nloc+6]);
        }
        break;
      default:
        cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
      return -1;
      }
      break;
    default:
      cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
      return -1;
    }
    
    return 0;
  }

  int partition(const vector<int> &ENList, const int& dim, int nloc, int nnodes, int npartitions, int partition_method, vector<int> &decomp){  
    // Build graph
    vector< set<int> > graph;
    int ret = FormGraph(ENList, dim, nloc, nnodes, graph);
    if(ret != 0){
      return ret;
    }
    
    // Compress graph    
    vector<idxtype> xadj(nnodes+1), adjncy;
    int pos=0;
    xadj[0]=1;
    for(int i=0;i<nnodes;i++){
      for(set<int>::iterator jt=graph[i].begin();jt!=graph[i].end();jt++){
        adjncy.push_back(*jt);
        pos++;
      }
      xadj[i+1] = pos+1;
      
      // Free memory as we go
      graph[i].clear();
    }
    graph.clear();
    
    // Partition graph
    decomp.resize(nnodes);
    int wgtflag=0, numflag=1, options[] = {0}, edgecut=0;
    
    if(partition_method){
      METIS_PartGraphKway(&nnodes, &(xadj[0]), &(adjncy[0]), NULL, NULL, &wgtflag, 
                          &numflag, &npartitions, options, &edgecut, &(decomp[0]));
    }else{
      METIS_PartGraphRecursive(&nnodes, &(xadj[0]), &(adjncy[0]), NULL, NULL, &wgtflag, 
                               &numflag, &npartitions, options, &edgecut, &(decomp[0]));
    }
    
    // number from zero
    for(int i=0;i<nnodes;i++)
      decomp[i]--;

    return edgecut;
  }  
  
  int partition(const vector<int> &ENList, int nloc, int nnodes, int npartitions, int partition_method, vector<int> &decomp){
    return partition(ENList, 3, nloc, nnodes, npartitions, partition_method, decomp);
  }
  
  int partition(const vector<int> &ENList, const vector<int> &surface_nids, const int& dim, int nloc, int nnodes, int npartitions, int partition_method, vector<int> &decomp){
    int num_elems = ENList.size()/nloc;
    
    set<int> surface_nodes;
    for(vector<int>::const_iterator it=surface_nids.begin(); it!=surface_nids.end(); ++it){
      surface_nodes.insert(*it);
    }

    // Build graph
    map<int, set<int> > graph;
    // vector< set<int> > graph(nnodes);
    switch (dim){
    case 3:
      switch (nloc){
      case 4:
        for(int i=0;i<num_elems;i++){
          for(int j=0;j<nloc;j++){
            if(surface_nodes.find(ENList[i*nloc+j])!=surface_nodes.end())
              for(int k=0;k<nloc;k++){
                if((j!=k)&&(surface_nodes.find(ENList[i*nloc+k])!=surface_nodes.end()))
                  graph[ENList[i*nloc+j]].insert(ENList[i*nloc+k]);
              }
          }
        }
        break;
      case 8:
        cerr<<"ERROR: Extrude support not implemented for hex's\n";
        exit(-1);
      default:
        cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
        return -1;
      }
      break;
    default:
      cerr<<"ERROR: element type not recognised - dim = "<<dim<<", nloc = "<<nloc<<endl; 
      return -1;
    }

    int snnodes=surface_nodes.size(); 
    
    // Compress graph
    vector<idxtype> xadj, adjncy;
    int pos=0;
    xadj.push_back(1);
    for(int i=0;i<snnodes;i++){
      assert(graph.find(i+1)!=graph.end());
      for(set<int>::iterator jt=graph[i+1].begin();jt!=graph[i+1].end();jt++){
        adjncy.push_back(*jt);
        pos++;
      }
      xadj.push_back(pos+1);
    }
    graph.clear();
    assert(xadj.size()==(snnodes+1));
    assert((*xadj.rbegin() - 1)==(adjncy.size()));

    // Partition graph
    vector<int> sdecomp(snnodes);
    int wgtflag=0, numflag=1, options[]={0}, edgecut=0;
    idxtype *vwgt=NULL, *adjwgt=NULL;
    if(partition_method){
      METIS_PartGraphKway(&snnodes, &(xadj[0]), &(adjncy[0]), vwgt, adjwgt, &wgtflag,
                          &numflag, &npartitions, options, &edgecut, &(sdecomp[0]));
    }else{
      METIS_PartGraphRecursive(&snnodes, &(xadj[0]), &(adjncy[0]), vwgt, adjwgt, &wgtflag,
                               &numflag, &npartitions, options, &edgecut, &(sdecomp[0]));
    }
    
    // -
    decomp.resize(nnodes);
    for(int i=0;i<nnodes;i++){
      decomp[i] = sdecomp[surface_nids[i]-1]-1;
    }
    
    return edgecut;
  }  
  
  int partition(const vector<int> &ENList, const vector<int> &surface_nids, int nloc, int nnodes, int npartitions, int partition_method, vector<int> &decomp){
                
    return partition(ENList, surface_nids, 3, nloc, nnodes, npartitions, partition_method, decomp);
  }
}
