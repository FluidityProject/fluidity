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
#include <mpi.h>
#include <string.h>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <set>
#include <assert.h>


#include "Mesh.h"
#include "samtypes.h"
#include "c++debug.h"

using namespace std;

void Mesh::export_fluidity(int int_memory[],         int int_memory_len,  int &IPT,
			   samfloat_t real_memory[], int real_memory_len, int &RPT,
			   const int MXNODS,
			   /* Based pointers into memory*/
			   int& NewPNodes, int& NewNNodes, int& NewNElems, int& NewNSElems, 
			   int& NewSCATER, int NewATOREC[],
			   int& NewGATHER, int NewATOSEN[],
			   int& NewENLIST, int& NewSNLIST, int& NewSURFID,
			   // Pressure stuff
			   int& NewPPressureNodes, int& NewNPressureNodes,
			   int& NewpSCATER,        int NewpATOREC[],
			   int& NewpGATHER,        int NewpATOSEN[],
			   int& NewpENLIST,
			   int& NewNODX, int& NewNODY, int& NewNODZ, samfloat_t Metric[]){ 
  ECHO("void Mesh::export_fluidity( ... )");
  CHECK(int_memory_len);
  CHECK(IPT);

  CHECK(real_memory_len);
  CHECK(RPT);

  //
  // Start exporting meshes
  //
  //  node_list.refresh_unn2gnn();

  NewPNodes = node_list.psize();
  NewNNodes = node_list.size();
  
  CHECK(NewPNodes);
  CHECK(NewNNodes);

  do_element_headcount();
  NewNElems = num_elements("volume");
  NewNSElems= num_elements("surface");

  // Start base pointer into the free integer memory.
  int bptr = IPT-1;
  
  {// Compress halo info..
    NewATOREC[0] = 1;
    for(int i=0; i<NProcs; i++){
      NewATOREC[i+1] =  NewATOREC[i] + halo_nodes[i].size();
    }
    
    NewSCATER = bptr + 1; // +1's is for bloody fortran
    for(int i=0; i<NProcs; i++){
      ECHO("Processor - " << i);
      for(set<unsigned>::iterator n = halo_nodes[i].begin(); n != halo_nodes[i].end(); ++n){
	int gnn = unn2gnn(*n);
	
	assert(gnn <  NewNNodes);
	assert(gnn >= NewPNodes);

	ECHO("SCATER = " << *n << " -> " << gnn + 1);
	
	int_memory[bptr++] = gnn + 1; // +1 is for the fortran offset.
      }
    }
  }
  
  {// Compress shared info...
    NewATOSEN[0] = 1;
    for(int i=0; i<NProcs; i++){
      NewATOSEN[i+1] = NewATOSEN[i] + shared_nodes[i].size();
    }
    
    NewGATHER = bptr + 1;
    for(int i = 0; i<NProcs; i++){
      ECHO("Processor - " << i);
      for(set<unsigned>::iterator n = shared_nodes[i].begin(); n != shared_nodes[i].end(); ++n){
	int gnn = unn2gnn(*n);
	
	assert(gnn <  NewPNodes);
	assert(gnn >= 0);

	ECHO("GATHER = " << *n << " -> " << gnn + 1);
	
	int_memory[bptr++] = gnn + 1; // +1 is for the fortran offset.
      }
    }
  } // end compressing shared
  
  { // Compress volume element-node lists
    NewENLIST = bptr + 1;
    for(deque<Element>::iterator et = element_list.begin(); et != element_list.end(); ++et){
      
      unsigned char type = (*et).get_flags();
      if( type & ELM_VOLUME ){
	vector<unn_t> enl( (*et).get_enlist() );
	for(vector<unn_t>::iterator it = enl.begin(); it != enl.end(); ++it){
	  int gnn = unn2gnn(*it);
	  
	  assert(bptr < int_memory_len );
	  assert(NewNNodes > gnn);
	  int_memory[bptr++] = gnn + 1;
	}
      }
    
    }
  }

  { 
    // Compress surface element-node list and write the surface id's.    
    vector<int> surfid( NewNSElems );
    unsigned pos=0;
    
    NewSNLIST = bptr + 1;
    for(deque<Element>::iterator et = element_list.begin(); et != element_list.end(); ++et){
      
      unsigned char type = (*et).get_flags();
      if( type & ELM_SURFACE ){
	
	assert(pos<NewNSElems);
	const vector<int>& ifields = (*et).get_ifields();
	surfid[pos++] = ifields[0];
	
	vector<unn_t> enl( (*et).get_enlist() );
	
	for(vector<unn_t>::iterator it = enl.begin(); it != enl.end(); ++it){
	  int gnn = unn2gnn(*it);
	  
	  assert(bptr < int_memory_len);
	  assert(NewNNodes > gnn);
	  int_memory[bptr++] = gnn + 1;
	}
      }
      
    }
    NewSURFID = bptr + 1;
    assert(NewNSElems == pos);
    for(size_t i=0; i<(size_t)NewNSElems; i++)
      int_memory[bptr++] = surfid[i];
    
  }

  // Export pressure nodes
  NewNPressureNodes = MFnode_list.size();
  NewPPressureNodes = MFnode_list.psize();
  CHECK(NewPPressureNodes);
  CHECK(NewNPressureNodes);
  
  if( mixed_formulation() ){
    ECHO("Doing mixed-formulation exporting thing.");
    
    {// Compress halo info..
      NewpATOREC[0]=1;
      for(int i = 0; i< NProcs; i++){
	NewpATOREC[i+1] = NewpATOREC[i] + halo_pnodes[i].size();
      }
      
      NewpSCATER = bptr + 1;
      for(int i=0; i<NProcs; i++){
	ECHO("Processor - " << i);
	for(set<unsigned>::iterator n = halo_pnodes[i].begin(); n != halo_pnodes[i].end(); ++n){
	  int gnn = MFunn2gnn(*n);
	  
	  assert(gnn <  NewNPressureNodes);
	  assert(gnn >= NewPPressureNodes);
	  ECHO("pSCATER = " << *n << " -> " << gnn + 1);
  
	  int_memory[bptr++] = gnn + 1; // +1 is for the fortran offset.
	}
      }
    } // finished compressing halo info.
    
    {// Compress shared info...
      NewpATOSEN[0]=1;
      for(int i=0; i<NProcs; i++){
	NewpATOSEN[i+1] = NewpATOSEN[i] + shared_pnodes[i].size();
      }
      
      NewpGATHER = bptr + 1;
      for(int i = 0; i<NProcs; i++){
	ECHO("Processor - " << i);
	for(set<unsigned>::iterator n = shared_pnodes[i].begin(); n != shared_pnodes[i].end(); ++n){
	  int gnn = MFunn2gnn(*n);
	  
	  assert(gnn <  NewPPressureNodes);
	  assert(gnn >= 0);
	  ECHO("pGATHER = " << *n << " -> " << gnn + 1);
	  
	  int_memory[bptr++] = gnn + 1; // +1 is for the fortran offset.
	}
      }
    } // end compressing shared
 
    ECHO("Writting pressure elements...");
    
    NewpENLIST = bptr + 1;
    for(ElementVector<Element>::iterator et = element_list.begin(); et != element_list.end(); ++et){
      if((*et).get_flags() & ELM_SURFACE)
	continue;
      
      const vector<unn_t>& enl = (*et).get_MFenlist();
      for(vector<unn_t>::const_iterator it = enl.begin(); it != enl.end(); ++it){
	int gnn = MFunn2gnn(*it);

	assert(bptr<int_memory_len);
	assert(NewNPressureNodes>gnn);
	int_memory[bptr++] = gnn + 1;
      }
    }
    ECHO("...done.");

  }

  IPT = bptr+1;

  // Start base pointer into the free real memory.
  { // Write coordinates, field data and metric.
    assert((3*MXNODS)<real_memory_len);
    bptr = RPT-1;
    RPT += MXNODS;
    NewNODX = bptr + 1;
    for(deque<Node>::iterator in=node_list.begin(); in != node_list.end(); ++in)
      real_memory[bptr++] = (*in).get_x();
    
    bptr = RPT-1;
    RPT += MXNODS;
    NewNODY = bptr + 1;  
    for(deque<Node>::iterator in=node_list.begin(); in != node_list.end(); ++in)
      real_memory[bptr++] = (*in).get_y();
    
    bptr = RPT-1;
    RPT += MXNODS;
    NewNODZ = bptr + 1;
    for(deque<Node>::iterator in=node_list.begin(); in != node_list.end(); ++in)
      real_memory[bptr++] = (*in).get_z();
    
    bptr = RPT-1;
    unsigned nfields = node_list[0].get_fields().size();
    RPT += (MXNODS*nfields);
    
    assert(RPT<real_memory_len);
    
    for(deque<Node>::iterator in=node_list.begin(); in != node_list.end(); ++in){
      CHECK(*in);
      
      const vector<samfloat_t>& flds = (*in).get_fields();
      assert(flds.size() == nfields);
      
      for(unsigned f=0; f<nfields; f++)
	real_memory[bptr + f*MXNODS] = flds[f];
      bptr++;
    }
    ECHO(RPT<<" should be the same as PYSFU2");

    int n=0;
    for(deque<Node>::iterator in=node_list.begin(); in != node_list.end(); in++,n++)
      memcpy(Metric+n*9, &(in->get_metric()[0]), 9*sizeof(samfloat_t));
    
    RPT = bptr+1;
  }

}

void Mesh::export_halo(int* colgat, int* atosen, int* scater, int* atorec, const int* ncolga, const int* nscate, const int* nprocs)
{
  // Input check
  assert(*ncolga == get_ncolga());
  assert(*nscate == get_nscate());
  assert(*nprocs == this->NProcs);
  
  atosen[0]= 0;
  int index = 0;
  for(int i = 0;i < this->NProcs;i++)
  {
    assert(i < *nprocs);
    atosen[i + 1] = atosen[i] + shared_nodes[i].size();
    for(set<unsigned>::iterator n = shared_nodes[i].begin(); n != shared_nodes[i].end(); ++n)
    {
      int gnn = unn2gnn(*n);	
      assert(index < *ncolga);
	    colgat[index++] = gnn + 1;
    }
  }

  atorec[0] = 0;
  index = 0;
  for(int i = 0;i < this->NProcs;i++)
  {
    assert(i < *nprocs);
    atorec[i + 1] = atorec[i] + halo_nodes[i].size();
    for(set<unsigned>::iterator n = halo_nodes[i].begin(); n != halo_nodes[i].end(); ++n)
    {
      int gnn = unn2gnn(*n);	
      assert(index < *nscate);
	    scater[index++] = gnn + 1;
    }
  }
  
  return;
}


void Mesh::export_phalo(int* pcolgat, int* patosen, int* pscater, int* patorec, const int* pncolga, const int* pnscate, const int* nprocs)
{
  // Input check
  assert(*pncolga == get_pncolga());
  assert(*pnscate == get_pnscate());
  assert(*nprocs == this->NProcs);
  
  patosen[0]= 0;
  int index = 0;
  for(int i = 0;i < this->NProcs;i++)
  {
    assert(i < *nprocs);
    patosen[i + 1] = patosen[i] + shared_pnodes[i].size();
    for(set<unsigned>::iterator n = shared_pnodes[i].begin(); n != shared_pnodes[i].end(); ++n)
    {
      int gnn = unn2gnn(*n);	
      assert(index < *pncolga);
	    pcolgat[index++] = gnn + 1;
    }
  }

  patorec[0] = 0;
  index = 0;
  for(int i = 0;i < this->NProcs;i++)
  {
    assert(i < *nprocs);
    patorec[i + 1] = patorec[i] + halo_pnodes[i].size();
    for(set<unsigned>::iterator n = halo_pnodes[i].begin(); n != halo_pnodes[i].end(); ++n)
    {
      int gnn = unn2gnn(*n);	
      assert(index < *pnscate);
	    pscater[index++] = gnn + 1;
    }
  }
  
  return;
}
