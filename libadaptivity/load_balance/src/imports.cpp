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
// Function to convert the mesh from FLUIDITY's matrix format
// to horizons internal format.

#include "samtypes.h"
#include "Mesh.h"
#include "c++debug.h"

#include <string.h>

#include <cassert>
#include <cmath>
#include <deque>
#include <set>
#include <vector>

using namespace std;

// Limit the scope of these two sets to this file
static set<int> deleted_vol_elems;
static set<int> deleted_sur_elems;

void Mesh::import_fluidity(const int dim, const int NNOD, const int NELM, const int NSELM, const int stride,
			   const int field_len,
			   const samfloat_t X[],
			   const samfloat_t Y[],
			   const samfloat_t Z[],
			   const samfloat_t metric[],
			   const samfloat_t field[],
			   const int ENLIST[], const int nloc,
			   const int SNLIST[], const int SURFID[], const int snloc,
			   const int ATOREC[], const int SCATER[],
			   const int ATOSEN[], const int GATHER[]){
  
  deleted_vol_elems.clear();
  deleted_sur_elems.clear();

  dimension = dim;

  // Map nodal information to internal format.
  ECHO("Reading nodes into internal format... ");
  size_t flen = field_len;
  for(size_t n=0; n<(unsigned)NNOD; n++){
    Node node;
    
    node.set_gnn( n );
    node.set_current_owner( MyRank );
    
    // Nodewise field values.
    for(size_t f = 0; f<flen; f++)
      node.append_field(field + stride*f + n, 1);
    
    // Position.
    switch(dim){
      case 3:
        node.set_coord(X[n], Y[n], Z[n]);
	break;
      case 2:
        node.set_coord(X[n], Y[n]);
	break;
      default:
        ERROR("Invalid dimension");
    }
    
    // Nodewise metric used for adaptivity. This is handled seperately
    // from fields as it can be used to control the
    // graph-partitioning.
    node.set_metric(metric + n*dim*dim, dim*dim);
    
    CHECK(node);
    add_node(node);
  }
  ECHO("Read.");

  //
  // Correct owner information for halo nodes.
  // (j = ATOREC[i]-1) is a base pointer to a list 
  // of nodes to be received from partition i 
  // SCATTER[j] - 1 is the node in question.
  //
  vector< vector<unsigned> > halo( NProcs ); // This is handy later
  for(unsigned p=0; p<(unsigned)NProcs; p++){
    
    for(int j=(ATOREC[p]-1);j<(ATOREC[p+1]-1);j++){
      unsigned haloNode = SCATER[j] - 1;

      assert(j<ATOREC[NProcs]-1);          // Crap out!
      assert(haloNode<(unsigned)NNOD);     // Crap out!
      
      node_list[ haloNode ].set_current_owner( p );
      node_list[ haloNode ].set_flags( NODE_HALO ) ;
      halo[p].push_back( haloNode );
    }
  }
  do_node_headcount();
  
  // Establish UNNing for native nodes
  {
    // Set up universal node numbering offset.
    find_nodes_per_rank( num_nodes("native") );
    
    unn_offsets[0] = 0;
    for(int i = 0; i<NProcs; i++)
      unn_offsets[i+1] = unn_offsets[i] + nodes_per_rank[i];
    
    unsigned cnt=0;
    for(unsigned i=0; i<(unsigned)NNOD; i++){
      if(node_list[i].get_current_owner() == MyRank){
	node_list[i].set_unn( unn_offsets[MyRank]+cnt );
	cnt++;
      }
    }
  }
  
  for(deque<Node>::iterator it=node_list.begin(); it!=node_list.end(); ++it){
    CHECK( (*it).get_gnn() );
    CHECK( (*it).get_current_owner() );
  }
  
  //
  // Record the nodes that are shared.
  // (j = ATOSEN[i]-1) is a base pointer to a list 
  // of nodes to be sent to partition i 
  // GATHER[j] - 1 is the node in question.
  //
  ECHO( "Recording shared nodes...." );

  // shared_unns is created to maintain the order of the shared
  // nodes so that the respective halo_node unn's can be correctly 
  // updated.
  vector< vector<unsigned> > shared_unns(NProcs);
  for(unsigned p=0; p<(unsigned)NProcs; p++){
    for(int j=(ATOSEN[p]-1);j<(ATOSEN[p+1]-1);j++){
      unsigned haloNode = GATHER[j] - 1;
      unsigned unn = node_list[haloNode].get_unn();
      
      // Store such that the halo is always in order of increasing unn.
      shared_nodes[p].insert(unn);
      shared_unns[p].push_back( unn );
    }
  }
  ECHO( "Recorded." );

  ECHO("At this point the node data should be converted to the internal format.");
  ECHO("Now it is time to correct the UNN's for the halo nodes.");
  
  //
  // Correct the UNN of the halo nodes (aka. SCATER nodes)
  //
  {
    vector<unsigned> halo_cnt( NProcs );
    for(int i = 0; i<NProcs; i++)
      halo_cnt[i] = halo[i].size();
    
    vector< vector<unsigned> > halo_unns;
    halo_update(shared_unns, halo_cnt, halo_unns);

    // Write in the UNN's and halo_nodes map
    for(unsigned i=0; i<(unsigned)NProcs; i++){

      assert( halo[i].size() == halo_unns[i].size() );
      
      const unsigned len = halo_unns[i].size();
            
      for(unsigned j = 0; j<len; j++){
	unsigned haloNode = halo[i][j];

	ECHO( i << ":" << haloNode << " --> " << halo_unns[i][j] );

	node_list[ haloNode ].set_unn( halo_unns[i][j] );
	halo_nodes[i].insert(halo_unns[i][j]);
      }
    }
    
  }

  ECHO( "gnn" <<"\t" << "unn" << "\t" << "owner" << "\t" << "unn2gnn");
  for(deque<Node>::iterator it=node_list.begin(); it!=node_list.end(); ++it){
    CHECK( (*it).get_gnn() );
    CHECK( (*it).get_current_owner() );
    CHECK( unn2gnn( (*it).get_unn() ) );
  }

  ECHO("UNN's have been established. ");
  ECHO("Lets check UNN's. Printing the unn's of the nodes in gather");
  
  for(unsigned i=0; i<(unsigned)NProcs; i++){
    ECHO("Rank " << i << ":\t");
    
    ECHO("SCATER");
    for(int j = (ATOREC[i]-1); j< (ATOREC[i+1]-1); j++)
      ECHO("\t\t(gnn="<<SCATER[j]-1 <<") "<< node_list[SCATER[j]-1].get_unn());
    ECHO("GATHER");
    for(int j = (ATOSEN[i]-1); j< (ATOSEN[i+1]-1); j++)
      ECHO("\t\t(gnn="<< GATHER[j]-1 <<")"<< node_list[GATHER[j]-1].get_unn());
  }
  
  unsigned eid=0;
  { // convert volume element information
    ECHO( "Write volume elements to internal format. " );
    CHECK(NELM);
    CHECK(nloc);

    Element elem;
    for(int e=0; e<NELM; e++){
      CHECK(e);
      elem.set_flags(ELM_DEFAULT|ELM_VOLUME);
      bool wanted = false; // Maybe it's not even this domain!
      
      { // write element node list with UNN's
	int cnt = nloc;
	int pos = e*nloc;                                // Note: no -1 here
	vector<unn_t> nl(cnt);
	for(int n=0; n<cnt; n++){
	  int nodeNum = ENLIST[pos + n] - 1;
	  assert(nodeNum>=0   ); // or crap out
	  assert(nodeNum<NNOD ); // or crap out
	  nl[ n ] = node_list[nodeNum].get_unn();
          if( node_list[nodeNum].get_current_owner() == MyRank )
            wanted = true;
	}
	elem.set_enlist( nl );
      }
      
      if( !wanted ){
	deleted_vol_elems.insert(e);
	continue;
      }

      // USING_DISCONTINUOUS
      // Field values??

      elem.set_eid(eid++);

      ECHO("Adding volume element:");
      CHECK(elem);
      
      add_element( elem );
    }
  }
  
  { // convert surface element information
    ECHO( "Write surface elements to internal format. " );
    
    Element elem;
    for(int e=0; e<NSELM; e++){
      elem.set_flags(ELM_DEFAULT|ELM_SURFACE);
      bool wanted = false; // Maybe it's on the halo interface

      { // write element node list with UNN's
	int cnt = snloc;
	int pos = e*snloc;                                // Note: no -1 here
	vector<unn_t> nl(cnt);
	for(int n=0; n<cnt; n++){
	  int nodeNum = SNLIST[pos + n] - 1;
	  assert(nodeNum>=0   ); // or crap out
	  assert(nodeNum<NNOD ); // or crap out
	  nl[ n ] = node_list[nodeNum].get_unn();
          if( node_list[nodeNum].get_current_owner() == MyRank )
            wanted = true;
	}
        if( !wanted ){
	  deleted_sur_elems.insert(e);
          continue;
	}

	elem.set_enlist( nl );
      }
      
      elem.set_eid(eid++);

      // add surface id's
      vector<int> ifield(1);
      ifield[0] = SURFID[e];
      elem.set_ifields( ifield );

      // USING_DISCONTINUOUS
      // Field values??
      ECHO( "Adding element:");
      CHECK( elem );

      add_element( elem );
    }
  do_element_headcount();
  }
 

  ECHO( "Finished translating fluidity stuff to internal format.");

  return;
}



void Mesh::import_pressure(const int NNOD, const int NELM,
			   const samfloat_t pressure[], const int nfields,
			   const int ENLIST[], const int ENLBAS[],
			   const int ATOREC[], const int SCATER[],
			   const int ATOSEN[], const int GATHER[]){
  // Make space for the pressure node list
  MFnode_list.clear();
  
  if(NNOD==0){
    
    // So the pressure mesh is the same as the velocity mesh so we can
    // just recreate it after we're finished.
    return;
    
  }else{
    
    // Set up UNNing for pressure nodes. As the number of pressure
    // nodes is always less or equal to the number of velocity nodes, and
    // the pressure node unn need not be contineous, the same offsets
    // will be used.    
    unsigned offset = 0;
    for(int i = 0; i<MyRank; i++)
      offset += nodes_per_rank[i];
    
    ECHO( "Reading pressure nodes into internal format... " );

    { // Map nodal information to internal format.
      for(unsigned n=0; n<(unsigned)NNOD; n++){
	PressureNode node;
	
	node.set_gnn( n );
	node.set_unn( offset + n );
	
	if(nfields>0){
	  vector<samfloat_t> flds(nfields);
	  memcpy(&(flds[0]), pressure+n*nfields, nfields*sizeof(samfloat_t));
	  node.set_pressure( flds );
	}
	
	MFnode_list.push_back( node );
      }
    }
    ECHO( "Read." );
    
    // shared_unns is created to maintain the order of the shared
    // nodes so that the respective halo_node unn's can be correctly 
    // updated.
    vector< vector<unsigned> > shared_unns(NProcs);
    for(unsigned p=0; p<(unsigned)NProcs; p++){
      for(int j=(ATOSEN[p]-1);j<(ATOSEN[p+1]-1);j++){
	
	unsigned haloNode = GATHER[j] - 1;
	unsigned unn = MFnode_list[haloNode].get_unn();
	
	// Store such that the halo is always in order of increasing unn.
	shared_unns[p].push_back( unn );
      }
    }
    
    //
    // Correct the UNN of the halo nodes (aka. SCATER nodes)
    //
    {    
      vector<unsigned> halo_cnt( NProcs );
      for(int i = 0; i<NProcs; i++)
	halo_cnt[i] = ATOREC[i+1] - ATOREC[i];
      
      vector< vector<unsigned> > halo_unns;
      halo_update(shared_unns, halo_cnt, halo_unns);
      
      // Write in the correct UNN's
      for(unsigned i=0; i<(unsigned)NProcs; i++){
	unsigned len = halo_unns[i].size();
	for(unsigned j = 0; j<len; j++){
	  unsigned haloNode = SCATER[ ATOREC[i]-1+j ];
	  MFnode_list[ haloNode ].set_unn( halo_unns[i][j] );
	}
      }
    }
    
    unsigned eid=0;    
    { // insert pressure node list into volume elements
      for(int e=0; e<NELM; e++){
	
	// Don't bother with elements that were deleted
	if(deleted_vol_elems.find(e) != deleted_vol_elems.end())
	  continue;
	
	int cnt = ENLBAS[e+1]-ENLBAS[e];
	int pos = ENLBAS[e];                                // Note: no -1 here
	vector<unn_t> nl(cnt);
	for(int n=0; n<cnt; n++){
	  int nodeNum = ENLIST[pos + n] - 1;
	  assert(nodeNum>=0   ); // or crap out
	  assert(nodeNum<NNOD ); // or crap out
	  nl[ n ] = MFnode_list[nodeNum].get_unn();
	}
	
	element_list[eid++].set_MFenlist( nl );
      }
    }
    
    //
    // Fluidity doesn't have surface pressure elements
    //
  }

  deleted_vol_elems.clear();
  deleted_sur_elems.clear();
}











