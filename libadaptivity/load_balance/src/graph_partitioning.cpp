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
#include "c++debug.h"
#include "Mesh.h"
#include "Graph.h"
#include "samtypes.h"

#include <vector>
#include <deque>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>

using namespace std;

std::vector<int> Mesh::decomp(const vector<int>& options){
  // options[0]:  x // Reserved for future use
  // options[1]:  1 // first partitioning
  //              2 // re-partitioning (diffusion)
  //              3 // re-partitioning (directed diffusion)
  //              4 // re-partitioning (clean)
  // options[2]:  1 // homogenious processing power
  //              2 // hetrogenious processing power
  // options[3]:  1 // no node weights
  //              2 // node weights based on the projected densitity of nodes
  // options[4]:  1 // no edge weights
  //              2 // edge weights designed for adaptivity

  // However the domain decomposition is going to be calculated, it is
  // going to be put into the vector noddom
  unsigned nnodes = num_nodes("total");
  vector<int> noddom( nnodes );

  ECHO("Building csr for ParMetis...."); 
  
  // Creating a node adjancy list in CSR format. This get initalised
  // using a node-node list generated using UNN's
  map<unsigned, set<unsigned> > edges = mknnlist(true);
  csr::Graph csr_graph( edges );  
  edges.clear();
  
  ECHO("...built.");
  
  // This may be required
  deque< set<unsigned> > nelist;

  { // Calculate any node weights
    CHECK( options[3] );
    switch( options[3] ){
    case 1:
      ECHO("No node weights");
      break;
    case 2:
      { // We wish to to limit the scope of a container
	ECHO("Optimize node weights for estimated node densitity");
	deque<float> elem_IdealDensity(__num_elements_total, -1.0);
	deque<float> elem_Volume(__num_elements_total, -1.0);
	if(nelist.empty())
	  nelist = mknelist();
	for(unsigned i=0; i<__num_elements_total; i++)
	  if(element_list[i].get_flags()&ELM_VOLUME){
	    elem_IdealDensity[i] = idealElementDensity(i);
	    elem_Volume[i]       = elementVolume(i);
	  }
	ECHO("Got element density field");

	ECHO("Calculate node weights");
	deque<float> weight(__num_nodes_total, 0.0);
	for(unsigned i=0;i<__num_nodes_total;i++){
	  for(set<unsigned>::const_iterator it=nelist[i].begin(); it!=nelist[i].end(); ++it){
	    if(elem_IdealDensity[*it] > 0.0){
	      weight[i] += elem_IdealDensity[*it];
	    }
	  }
	  weight[i]/=( 4*nelist[i].size() ); // lumped mass matrix approximation
	}
	ECHO("Finished calculating node weights");
	csr_graph.setVertexWeights(weight);
      }
      break;
    default:
      ERROR("Unknown option " << options[3] << " for node weights passed down into Sam. Giving up!");
      MPI::COMM_WORLD.Abort( MPI_ERR_OP );
      break;
    }
    ECHO("Got node weights.");    
  }
  
  
  { // Calculate any edge weights
    // Note: The really strange scoping came about because of some odd
    // compiler bug in gcc. Namely "crosses initialization of".    
    CHECK( options[4] );

    switch( options[4] ){
    case 1:
      ECHO("No edge weights");
      break;
    case 2:
      { // We wish to to limit the scope of a container
	ECHO("Optimize edge weights for adaptivity");
	deque<samfloat_t> elem_fxnls( __num_elements_total );
	if(nelist.empty())
	  nelist = mknelist();
	
	ECHO("Calculating all the element functionals...");
	for(unsigned i=0; i<__num_elements_total; i++)
	  elem_fxnls[i] = element_functional(i);
	ECHO("...done");
	
	{
	  ECHO("Calculating a unn2gnn mapping...");
	  map<unsigned, unsigned> unn2gnn;
	  unsigned pos = 0;
	  for(NodeVector<Node>::const_iterator it = node_list.begin(); it != node_list.end(); ++it){
	    unn2gnn[ (*it).get_unn() ] = pos++;
	  }
	  ECHO("...done");
	  
	  ECHO("Deriving edge weights...");
	  csr_graph.buildEdgeWeights_adapt(nelist, unn2gnn, elem_fxnls, 
					   functional_tolerence);
	  ECHO("...done");
	}
      }

      break;
    case 3: // Optimize edge weights for solver
      break;
    default:
      ERROR("Unknown option " << options[4] << " for edge weights passed down into Sam. Giving up!");
      MPI::COMM_WORLD.Abort( MPI_ERR_OP );
    }
    ECHO("Got edge weights.");
  }

  //
  // Check for heterogeneous computational power
  //
  vector<samfloat_t> flops_profile;
  if( options[2] == 2 ){ // Call the nurse

    assert(false);
#if 0
    ECHO("Calling NURSE.");

    // Make space to store times
    vector<samfloat_t> times(NUMBER_OF_BENCHMARK_TRIALS*NProcs);
    
    // Make space to store standard errors
    vector<samfloat_t> std_err(NProcs);
  
    // Make space to store means
    vector<samfloat_t> mean(NProcs);

    // Do trials
    for(unsigned trial=0; trial<NUMBER_OF_BENCHMARK_TRIALS; trial++){
      MPI_Barrier(MPI_COMM_WORLD);
      mpi_nbm(nbm_cacheater, times.begin() + trial*NProcs);
    }
    
    // Get mean
    for(unsigned p=0; p<(unsigned)NProcs; p++){
      samfloat_t sum = 0;
      for(unsigned i=0; i<NUMBER_OF_BENCHMARK_TRIALS; i++){
	sum += times[i*NProcs + p];
      }
      mean[p] = sum/NUMBER_OF_BENCHMARK_TRIALS;
    }

    // Get standard error
    for(unsigned p=0; p<(unsigned)NProcs; p++){
      samfloat_t sum = 0;
      for(unsigned i=0; i<NUMBER_OF_BENCHMARK_TRIALS; i++)
	sum += ( (times[i*NProcs + p] - mean[p])*(times[i*NProcs + p] - mean[p]) );
      std_err[p] = sqrt(sum/NUMBER_OF_BENCHMARK_TRIALS);
    }

    // Decide if it's heterogeneous or not.
    samfloat_t sum = 0.0;
    samfloat_t maxt = times[0];
    samfloat_t mint = times[0];
    for(unsigned p=0; p<(unsigned)NProcs; p++){      
      maxt = max(maxt, times[p]);
      mint = min(mint, times[p]);
      sum += times[p];
    }
    samfloat_t average = sum/NProcs;
    if( ((maxt-mint)/average) > HETERO_TOL ){
      heterogeneous_mflops = true;
      flops_profile.resize(NProcs);
      samfloat_t lnorm = 0.0;
      for(unsigned p=0; p<(unsigned)NProcs; p++){
	flops_profile[p] = 1.0 / times[p];
	lnorm += flops_profile[p];
      }
      for(unsigned p=0; p<(unsigned)NProcs; p++)
	flops_profile[p] /= lnorm;
    }
#endif
  } // Finished setting up heterogeneous support for processing power
  

  // Set up vtxdist.
  ECHO("size vtxdist = " << unn_offsets.size());
  for(int p=0; p<=NProcs; p++)
    ECHO("vtxdist["<<p<<"] = " <<  unn_offsets[p]);
  for(unsigned i=0; i<5; i++)
    ECHO("decomposition options["<<i<<"] = " << options[i]);
  
  
  // Chose your domain decomposition method
  switch( options[1] ){
  case 1: // clean repartitioning into nparts
    
    ECHO("Doing ParMETIS V3 PartKway .. ");
    csr_graph.Repart(unn_offsets, noddom, options[0]);
    ECHO("...done.");
    break;
    
  default:
    
    ECHO("Doing RepartRemap (ParMetis) .. ");
    csr_graph.RepartRemap(unn_offsets, noddom);
    ECHO("...done.");
    break;
  }
  
  ECHO("Finished ParMetis.");

  //
  // ParMetis has only written in the node distribution for the native nodes.
  // Write in the new distribution for the halo as well.
  //
  vector< vector<unsigned> > shared_noddom(NProcs);
  for(int i=0; i<NProcs; i++){
    int len = shared_nodes[i].size();
    if(len==0) 
      continue;

    shared_noddom[i].resize( len );
    unsigned cnt=0;
    for(set<unsigned>::iterator it = shared_nodes[i].begin(); it != shared_nodes[i].end(); ++it)
      shared_noddom[i][cnt++] = noddom[unn2gnn(*it)];
  }
  
  vector< vector<unsigned> > halo_noddom;
  halo_update(shared_noddom, halo_noddom);
  
  // Write in the mapping
  for(unsigned i=0; i<(unsigned)NProcs; i++){
    unsigned cnt=0;
    for(set<unsigned>::iterator it = halo_nodes[i].begin(); it != halo_nodes[i].end(); ++it){
      noddom[unn2gnn(*it)] = halo_noddom[i][cnt++];
    }
  }
  
  return noddom;
}
