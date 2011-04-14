/* Copyright (C) 2006 Imperial College London and others.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

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
//
// provides interface for ParMetis
//
#include <mpi.h>
#include <vector>
#include <set>
#include <deque>
#include <cstdlib>
#include <assert.h>
#include <string.h>
#include <iostream>

#include "Graph.h"
#include "samParMetis.h"
#include "c++debug.h"

using namespace std;
using namespace csr;

void Graph::RepartRemap(vector<int>& vtxdist, vector<int>& noddom){
  MPI_Comm comm_world = MPI_COMM_WORLD;
  int nparts = MPI::COMM_WORLD.Get_size();

  int *vwgt   = NULL;
  int *adjwgt = NULL;
  
  if(!nweight.empty())
    vwgt = &(nweight[0]);
  if(!eweight.empty())
    adjwgt = &(eweight[0]);
  
  int wgtflag = 0;
  if((vwgt==NULL)&&(adjwgt!=NULL)){
    wgtflag = 1;
  }else if((vwgt!=NULL)&&(adjwgt==NULL)){
    wgtflag = 2;
  }else if((vwgt!=NULL)&&(adjwgt!=NULL)){
    wgtflag = 3;
  }

  int ncon = 1;
  vector<float> tpwgts(nparts*ncon), ubvec(ncon);
  for(int i=0;i<nparts*ncon;i++)
    tpwgts[i] = 1.0/(float)(nparts);
  for (int i=0;i<ncon;i++)
    ubvec[i] = 1.05;
  
  float ipc_factor = 100000.0;
  int options[4];
  options[0] = 0; // not using options

  int edgecut;
  int numflag = 0;
  
  ParMETIS_V3_AdaptiveRepart(&(vtxdist[0]), &(bptr[0]), &(edges[0]), vwgt, NULL,
                             adjwgt, &wgtflag, &numflag, &ncon, &nparts, &(tpwgts[0]),
                             &(ubvec[0]), &ipc_factor, options, &edgecut, &(noddom[0]), &comm_world);

  return;
}

void Graph::Repart(vector<int>& vtxdist, vector<int>& noddom, int nparts){
  int edgecut;
  int options[4];
  int numflag = 0;
  
  int nprocs = MPI::COMM_WORLD.Get_size();
  int MyRank = MPI::COMM_WORLD.Get_rank();

  vector<int> nnodes(nprocs, 0);
  nnodes[MyRank] = noddom.size();
  
  MPI::COMM_WORLD.Allgather(&(nnodes[MyRank]), 1, MPI_INT, &(nnodes[0]), 1, MPI_INT);
  
  vector<int> ranks;
  for(int i=0;i<nprocs;i++){
    if(nnodes[i])
      ranks.push_back(i);
  }

  if(ranks.size()==0){
    cerr<<"ERROR: None of the domains have any nodes. There is apparently no problem"
        <<" so I'm out of here.\n";
    exit(-1);
  }

  MPI_Group world_group, sub_group; 
  MPI_Comm world_comm=MPI_COMM_WORLD, sub_comm; 
  MPI_Comm_group(world_comm, &world_group); 
  MPI_Group_incl(world_group, ranks.size(), &(ranks[0]), &sub_group);
  MPI_Comm_create(world_comm, sub_group, &sub_comm); 
  
  int wgtflag = 0;
  if(eweight.empty()==false){
    wgtflag = 1;
  }
  
  options[0] = 0;
  
  int ncon=0;

  if(nnodes[MyRank]){
    ParMETIS_V3_PartKway(&(vtxdist[0]), &(bptr[0]), &(edges[0]), NULL, &(eweight[0]),
                         &wgtflag, &numflag, &ncon, &nparts, NULL, NULL,
                         options, &edgecut, &(noddom[0]), &sub_comm);

    MPI_Comm_free(&sub_comm);
  }
}






































