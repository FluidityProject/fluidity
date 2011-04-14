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

#include "c++debug.h"
#include "packing.h"
#include "Mesh.h"
#include "comTools.h"

#include <cassert>
#include <string>
#include <vector>

using namespace std;

static int packing_ncalls = 0;

packing::packing(const MI5& intelligance_report, const Mesh& mesh){
  
  NProcs   = MPI::COMM_WORLD.Get_size();
  MyRank   = MPI::COMM_WORLD.Get_rank();
  __nnodes = intelligance_report.nnodes();
  __nelems = intelligance_report.nelems();
  
  SendRecvBuffer.resize( NProcs );
  
  size_t max_nodepack_size     = mesh.max_nodepack_size();
  size_t max_elempack_size     = mesh.max_elementpack_size();
  size_t max_pressurepack_size = mesh.max_pressurepack_size();

  size_t max_nodes2send = intelligance_report.nodes2send[0].size();
  size_t max_elems2send = intelligance_report.elems2send[0].size();

  size_t sum_nodes2send = 0;
  size_t sum_elems2send = 0;

  // allocate space for the send buffers.
  size_t nbytes;
  int _space_for_ints = MPI::INT.Pack_size(3, MPI::COMM_WORLD);

  for(unsigned p=0; p<NProcs; p++){

    max_nodes2send = max(max_nodes2send, intelligance_report.nodes2send[p].size());
    max_elems2send = max(max_elems2send, intelligance_report.elems2send[p].size());

    sum_nodes2send += intelligance_report.nodes2send[p].size();
    sum_elems2send += intelligance_report.elems2send[p].size();
    
    nbytes = _space_for_ints +
      (intelligance_report.nodes2send[p].size())*max_nodepack_size + 
      (intelligance_report.pnodes2send[p].size())*max_pressurepack_size + 
      (intelligance_report.elems2send[p].size())*max_elempack_size;

    // For good measure
    nbytes = (unsigned)(nbytes*1.1);
    
    ECHO("Requesting "<<nbytes<<" bytes for send buffer "<<p);
    
    SendRecvBuffer[p].resize(nbytes);
  }
  
//  cout<<"migrating"<<packing_ncalls<<", "
//      <<max_nodes2send<<" "<<max_elems2send<<" "
//      <<sum_nodes2send<<" "<<sum_elems2send<<endl;

  // --
  offsets.resize( NProcs );
  for(unsigned p=0; p<NProcs; p++) offsets[p] = 0;

  packing_ncalls++;
}

packing::~packing(){
  clear();
}

void packing::pack(const MI5& intelligance_report, Mesh& mesh){
  char *buff;
  
  ECHO("Starting to pack");

  //
  // Pack node data
  //
  for(unsigned p=0; p<NProcs; p++){
    int len = SendRecvBuffer[p].size();
    buff    = &(SendRecvBuffer[p][0]);
    int cnt = intelligance_report.nodes2send[p].size();

    ECHO("Packing "<<cnt<<" nodes for "<<p);
    
    MPI::INT.Pack(&cnt, 1, buff, len, offsets[p], MPI::COMM_WORLD);
  }
  pack_nodes(intelligance_report, mesh);
  
  //
  // Pack element data
  //
  for(unsigned p=0; p<NProcs; p++){
    int len = SendRecvBuffer[p].size();
    buff    = &(SendRecvBuffer[p][0]);
    int cnt = intelligance_report.elems2send[p].size();

    ECHO("Packing "<<cnt<<" elements for "<<p);
    
    MPI::INT.Pack(&cnt, 1, buff, len, offsets[p], MPI::COMM_WORLD);
  }
  pack_elems(intelligance_report, mesh);
  
  //
  // Pack pressure-node data
  //
  for(unsigned p=0; p<NProcs; p++){
    int len = SendRecvBuffer[p].size();
    buff    = &(SendRecvBuffer[p][0]);
    int cnt = intelligance_report.pnodes2send[p].size();
    
    ECHO("Pressure "<<cnt<<" nodes for "<<p<<" buff len="<<len<<", offsets[p]="<<offsets[p]);
    
    MPI::INT.Pack(&cnt, 1, buff, len, offsets[p], MPI::COMM_WORLD);
  }
  pack_pnodes(intelligance_report, mesh);
  
  // 
  // We might possibly resize the send buffers here but I still don't
  // know if this is a good idea or not.
  //
#ifdef MINIMIZE_MEMORY
  for(unsigned p=0; p<NProcs; p++){
    SendRecvBuffer[p].resize( offsets[p] );
  }
#endif

  ECHO("Packing Complete");
}

void packing::pack_nodes(const MI5& intelligance_report, Mesh& mesh){
  char *buffer;
  
  for(unsigned p=0; p<NProcs; p++){
    ECHO("Node pack for "<< p);

    unsigned cnt = intelligance_report.nodes2send[p].size();
    for(unsigned n=0; n<cnt; n++){
      
      int node = intelligance_report.nodes2send[p][n];
      int bsize = SendRecvBuffer[p].size();      
      buffer    = &(SendRecvBuffer[p][0]);

      Node __node__ = mesh.get_node(node);
      
      __node__.pack(buffer, bsize, offsets[p]);
    }
  }
}

void packing::pack_pnodes(const MI5& intelligance_report, Mesh& mesh){
  char *buffer;
  
  for(unsigned p=0; p<NProcs; p++){
    ECHO("Pressure-node pack for "<< p);
    
    unsigned cnt = intelligance_report.pnodes2send[p].size();
    for(unsigned n=0; n<cnt; n++){
      
      int node = intelligance_report.pnodes2send[p][n];
      int bsize = SendRecvBuffer[p].size();      
      buffer    = &(SendRecvBuffer[p][0]);
      
      const PressureNode& __node__ = mesh.get_pnode(node);
      
      __node__.pack(buffer, bsize, offsets[p]);
    }
  }
}

void packing::pack_elems(const MI5& intelligance_report, Mesh& mesh){
  char *buffer;

  for(unsigned p=0; p<NProcs; p++){
    
    // Set up the buffer.
    int bsize = SendRecvBuffer[p].size();      
    buffer    = &(SendRecvBuffer[p][0]);
    
    for(unsigned e=0; e<intelligance_report.elems2send[p].size(); e++){
      
      int elem = intelligance_report.elems2send[p][e];      

      Element __elem__ = mesh.get_element(elem);
      __elem__.pack(buffer, bsize, offsets[p]);

    }
  }
}

void packing::unpack(MI5& intelligance_report, Mesh& mesh){

  // Reset offsets
  for(unsigned p=0; p<NProcs; p++)
    offsets[p] = 0;
  
  // Unpack node data.
  unpack_nodes(intelligance_report, mesh);
  
  // Unpack element data
  unpack_elems(intelligance_report, mesh);

  // Unpack pressure-node data.
  unpack_pnodes(intelligance_report, mesh);

}

void packing::unpack_nodes( MI5& intelligance_report, Mesh& mesh){
  char *buffer;
  
  ECHO("Unpacking nodes!");
  
  for(unsigned p=0; p<NProcs; p++){
    if( p == MyRank) 
      continue;

    int ncnt;
    int nbytes = SendRecvBuffer[p].size();
    buffer     = &(SendRecvBuffer[p][0]);

    // How many new nodes.
    MPI::INT.Unpack(buffer, nbytes, &ncnt, 1, offsets[p], MPI::COMM_WORLD);
    
    ECHO("Unpacking " << ncnt << " nodes from " << p);
    
    // --
    for(unsigned n = 0; n<(unsigned)ncnt; n++){
      Node node;
      
      // Unpack node...
      node.unpack(buffer, nbytes, offsets[p]);
     
      ECHO("Unpacking " << ncnt);
      CHECK(node);

      // update new_owned_nodes and new_halo_nodes
      if( node.get_future_owner() == MyRank ){  // New owned...
	intelligance_report.new_owned_nodes.push_back( mesh.num_nodes("total")+n );
	mesh.add_node( node );
      }else{ // new halo
	intelligance_report.new_halo_nodes[ node.get_future_owner() ][node.get_unn()] = node;
      }
    }
  }
}

void packing::unpack_elems(MI5& intelligance_report, Mesh& mesh){
  char *buffer;
  
  for(unsigned p=0; p<NProcs; p++){
    if(p==MyRank)
      continue;

    int ecnt;

    int nbytes = SendRecvBuffer[p].size();
    buffer     = &(SendRecvBuffer[p][0]);
    
    // How many new elements.
    MPI::INT.Unpack(buffer, nbytes, &ecnt, 1, offsets[p], MPI::COMM_WORLD);

    ECHO("Unpacking " << ecnt << " elements from " << p);

    for(int elem=0; elem<ecnt; elem++){
      Element element;
      
      // Unpack element...      
      element.unpack(buffer, nbytes, offsets[p]);

      ECHO("Unpacking " << ecnt);
      CHECK(element);
  
      mesh.add_element(element);
    }
  }
}

void packing::unpack_pnodes( MI5& intelligance_report, Mesh& mesh){
  ECHO("Unpacking pressure nodes...");
  
  for(unsigned p=0; p<NProcs; p++){
    if( p == MyRank) 
      continue;
    
    int ncnt;
    int nbytes = SendRecvBuffer[p].size();
    char *buffer     = &(SendRecvBuffer[p][0]);
    
    // How many new pressure nodes.
    MPI::INT.Unpack(buffer, nbytes, &ncnt, 1, offsets[p], MPI::COMM_WORLD);
    
    ECHO("Unpacking " << ncnt << " pressure nodes from " << p);
    
    // --
    for(unsigned n = 0; n<(unsigned)ncnt; n++){
      PressureNode node;
      
      // Unpack node...
      node.unpack(buffer, nbytes, offsets[p]);
      
      ECHO("Unpacking " << ncnt);
      CHECK(node);
      
      mesh.add_node( node );
    }
  }

  ECHO("...unpacked!");
}

void packing::send(){
  
  ECHO("Sending submeshes.");
  allSendRecv(SendRecvBuffer);
  ECHO("Received submeshes.");
  
  return;
}

void packing::clear(){
  offsets.clear();
  SendRecvBuffer.clear();
}





