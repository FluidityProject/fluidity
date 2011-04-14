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
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <set>
#include <cassert>
#include <cstdlib>

#include "Mesh.h"
#include "c++debug.h"

#include "ElementVector.h"
#include "NodeVector.h"
#include "functional_2d.h"

using namespace std;

// Mesh class constructer. Initalize some variables.
Mesh::Mesh(){
  __mixed_formulation     = false;
  __num_nodes_native      = 0;
  __num_nodes_alien       = 0;
  __num_nodes_total       = 0;
  __num_elements_total    = 0;
  __num_elements_volume   = 0;
  __num_elements_surface  = 0;

  dimension = -1;
  
  if(MPI::Is_initialized()){
    MyRank = MPI::COMM_WORLD.Get_rank();
    NProcs = MPI::COMM_WORLD.Get_size();
  }else{
    MyRank = 0;
    NProcs = 1;
  }
  
  shared_nodes.resize( NProcs );
  halo_nodes.resize( NProcs ); 
  
  nodes_per_rank.resize(NProcs);
  unn_offsets.resize(NProcs+1);
}

Mesh::~Mesh(){}

// returns true if we're using mixed formulation
bool Mesh::mixed_formulation() const{
  return( __mixed_formulation );
}               
bool Mesh::mixed_formulation(const bool in){
  __mixed_formulation = in;
  return( __mixed_formulation );
}

// Append a node to our list.
void Mesh::add_node(const Node& node){
  node_list.push_back(node);
}
void Mesh::add_node(const PressureNode& node){
  MFnode_list.push_back(node);
}

void Mesh::do_node_headcount(){
  __num_nodes_native = 0;
  __num_nodes_alien  = 0;

  for(NodeVector<Node>::iterator it=node_list.begin(); 
      it != node_list.end(); ++it){
    if( (*it).get_current_owner() == (unsigned)MyRank){
      __num_nodes_native++;
    }else{
      __num_nodes_alien++;
    }
  }
  
  __num_nodes_total = __num_nodes_native + __num_nodes_alien;
  assert( __num_nodes_total == node_list.size());
  
}
void Mesh::do_element_headcount(){
  __num_elements_surface = 0;
  __num_elements_volume  = 0;

  for(ElementVector<Element>::iterator it = element_list.begin(); 
      it != element_list.end(); ++it){
    const unsigned char flag = (*it).get_flags();

    if( flag & ELM_SURFACE){
      __num_elements_surface++;
    }else if( flag & ELM_VOLUME){
      __num_elements_volume++;
    }else{
      ERROR("crap! Unreconised element type");
      exit(-2); // crap out
    }
  }

  __num_elements_total = __num_elements_volume + __num_elements_surface;
  assert( __num_elements_total == element_list.size() );
}

unsigned Mesh::max_nodepack_size() const{
  unsigned maxsize = 0;
  
  for(NodeVector<Node>::const_iterator nt=node_list.begin(); nt!=node_list.end(); ++nt){
    unsigned size = (*nt).pack_size();
    maxsize = (size>maxsize)?size:maxsize;
  }
  
  return(maxsize);
}

unsigned Mesh::max_pressurepack_size() const{
  unsigned maxsize = 0;
  
  for(NodeVector<PressureNode>::const_iterator nt=MFnode_list.begin(); nt!=MFnode_list.end(); ++nt){
    unsigned size = (*nt).pack_size();
    maxsize = (size>maxsize)?size:maxsize;
  }
  
  return(maxsize);
}


void Mesh::add_element(Element& elem){
  element_list.push_back(elem);
}

unsigned int Mesh::max_elementpack_size() const{
  unsigned maxsize = 0;
  
  for(ElementVector<Element>::const_iterator et=element_list.begin(); et!=element_list.end(); ++et){
    unsigned size = (*et).pack_size();
    maxsize = (size>maxsize)?size:maxsize;
  }
  
  return(maxsize);
}

void Mesh::set_node_current_owner(const unsigned int n, const unsigned short o){
  node_list[n].set_current_owner(o);
}

unsigned Mesh::get_node_current_owner(const unsigned n){
  return node_list[n].get_current_owner();
}

void Mesh::set_node_future_owner(const unsigned int n, const unsigned short o){
  node_list[n].set_future_owner(o);
}

unsigned Mesh::get_node_future_owner(const unsigned int n){
  return node_list[n].get_future_owner();
}

void Mesh::set_node_flags(const unsigned int n, const unsigned char f){
  node_list[n].set_flags( f );
}

void Mesh::set_node_unn(const unsigned int n, const unn_t unn){
  node_list[n].set_unn( unn );
}

unn_t Mesh::get_node_unn(const unsigned int n){
  return node_list[n].get_unn();
}

std::vector<unn_t> Mesh::get_enlist(const unsigned int e){
  return element_list[e].get_enlist();
}

unsigned int Mesh::get_enlist_size(const unsigned int e){
  return element_list[e].get_size_enlist();
}

const Node& Mesh::get_node(const unsigned n) const{
  return node_list[n];
}
const PressureNode& Mesh::get_pnode(const unsigned n) const{
  const PressureNode& node = MFnode_list[n];
  return node;
}
  
Element& Mesh::get_element(const unsigned e){
  return element_list[e];
}

unsigned char Mesh::get_node_flags(const unsigned int n){
  return node_list[n].get_flags();
}
vector<samfloat_t> Mesh::get_node_fields(const unsigned int n){
  return node_list[n].get_fields();
}

gnn_t Mesh::unn2gnn(const unn_t unn){
  const Node& node = node_list.unn(unn);
  return node.get_gnn();
}
gnn_t Mesh::MFunn2gnn(const unn_t unn){
  const PressureNode& node = MFnode_list.unn(unn);
  return node.get_gnn();
}

unsigned Mesh::num_nodes(const string option){
  if(option == "native")
    return __num_nodes_native;
  if(option == "alien")
    return __num_nodes_alien;
  
  if(option != "total"){
    ERROR( "Error: Invalid option, " << option << 
	   ", passed to unsigned Mesh::num_nodes(const string option)" );
  }
  
  return __num_nodes_total;
  
}

unsigned Mesh::num_elements(const string& option){
  if(option == "surface"){
    return __num_elements_surface;
  }else if(option == "volume"){
    return __num_elements_volume;
  }else if(option == "total"){
    return __num_elements_total;
  }else{
    assert(1); // Crap out!
  }

  return 0; // just so that we get no warnings from the compiler
}

unsigned Mesh::num_nodes_shared(){
  unsigned int cnt=0;
  
  for(int i = 0; i<NProcs; i++)
    cnt += shared_nodes[i].size();
  
  return( cnt );
}
unsigned Mesh::num_nodes_shared(const unsigned proc){
  return shared_nodes[ proc ].size();
}

unsigned Mesh::num_nodes_halo(){
  unsigned int cnt=0;
  
  for(int i = 0; i<NProcs; i++)
    cnt += halo_nodes[i].size();
  
  return cnt;
}
unsigned Mesh::num_nodes_halo(const unsigned proc){
  return halo_nodes[proc].size();
}

void Mesh::find_nodes_per_rank(const unsigned lcnt){
  if( MPI::Is_initialized() ){
    nodes_per_rank.resize( NProcs );
    
    for(int i = 0; i<NProcs; i++)
      MPI::COMM_WORLD.Gather(&lcnt, 1, MPI::INT, &(nodes_per_rank[0]), 1, MPI::INT, i);
    
  }else{ // Not running MPI.
    nodes_per_rank.resize(1);
    nodes_per_rank[0] = lcnt;
  }
}

unsigned Mesh::num_nodes_in_rank(const unsigned rank){
  return nodes_per_rank[rank];
}

void Mesh::add2shared_nodes(const unsigned short p, const unn_t unn){
  shared_nodes[p].insert(unn);
}

void Mesh::add2halo_nodes(const unsigned short p, const unn_t unn){
  halo_nodes[p].insert(unn);
}

// find the set of domains which know this element.
set<int> Mesh::current_element_domain_set(const unsigned elem){
  set<int> doms;
  vector<unn_t> nodes = element_list[elem].get_enlist();
  
  for(vector<unn_t>::iterator it = nodes.begin(); it != nodes.end(); ++it){
    gnn_t gnn = unn2gnn(*it);

    assert( gnn < (gnn_t)num_nodes("total") );
    doms.insert( node_list[ gnn ].get_current_owner() );
  }
  
  return doms;
}
 
// find the set of domains which know this element.
set<int> Mesh::future_element_domain_set(const unsigned elem){
  set<int> doms;
  vector<unn_t> nodes = element_list[elem].get_enlist();
  for(vector<unn_t>::iterator it = nodes.begin(); it != nodes.end(); ++it){
     gnn_t gnn = unn2gnn(*it);

     assert( gnn < (gnn_t)num_nodes("total") );
     doms.insert( node_list[ gnn ].get_future_owner() );
  }

  return doms;
}
// current mapping of the nodes in element x to domains
map<unn_t, int> Mesh::current_enode2domain_map(const unsigned elem){
  map<unn_t, int> mapping;
  
  vector<unn_t> nodes = element_list[elem].get_enlist();
  for(vector<unn_t>::iterator it = nodes.begin(); it != nodes.end(); ++it){
    ECHO("unn? " << *it);
    mapping[ *it ] = node_list[ unn2gnn(*it) ].get_current_owner();
  }

  return mapping;
}
// future mapping of the nodes in element x to domains
map<unn_t, int> Mesh::future_enode2domain_map(const unsigned elem){
  map<unn_t, int> mapping;
  
  vector<unn_t> nodes = element_list[elem].get_enlist();
  for(vector<unn_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    mapping[ *it ] = node_list[ unn2gnn(*it) ].get_future_owner();
  
  return mapping;
}

// Updata all halo values stored in data_in
void Mesh::halo_update(const vector< vector<int> >& data_in, vector< vector<int> >& data_out){
  
  data_out.clear();
  data_out.resize(NProcs);

  // Calculate the number of items being sent per node
  unsigned stride=0;
  for(unsigned i=0;i<(unsigned)NProcs;i++){
    if(data_in[i].empty())
      continue;
    else{
      stride = data_in[i].size() / num_nodes_shared(i);
      break;
    }
  }
  
  // setup non-blocking receives
  vector<MPI::Request> RecvRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_halo(i);
    
    if((i==MyRank)||(cnt==0)){
      RecvRequest[i] =  MPI::REQUEST_NULL;
    }else{
      data_out[i].resize(cnt);
      RecvRequest[i] = MPI::COMM_WORLD.Irecv(&(data_out[i][0]), cnt, MPI::INT, i, 13);
    }
  }
  
  // setup non-blocking sends
  vector<MPI::Request> SendRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_shared(i);

    if((i==MyRank)||(cnt==0)){
      SendRequest[i] =  MPI::REQUEST_NULL;
    }else{    
      SendRequest[i] = MPI::COMM_WORLD.Isend(&(data_in[i][0]), cnt, MPI::INT, i, 13);
    }
  }
  
  MPI::Request::Waitall(NProcs, &(RecvRequest[0]));
  MPI::Request::Waitall(NProcs, &(SendRequest[0]));
  
}
// Updata all halo values stored in data_in
void Mesh::halo_update(const vector< vector<unsigned> >& data_in, vector< vector<unsigned> >& data_out){
  
  data_out.clear();
  data_out.resize(NProcs);

  // Calculate the number of items being sent per node
  unsigned stride=0;
  for(unsigned i=0; i<(unsigned)NProcs; i++){
    if(data_in[i].empty())
      continue;
    else{
      stride = data_in[i].size() / num_nodes_shared(i);
      break;
    }
  }
  
  // setup non-blocking receives
  vector<MPI::Request> RecvRequest(NProcs);
  for(int i=0;i<NProcs;i++){
    int cnt = stride*num_nodes_halo(i);

    if((i==MyRank)||(cnt==0)){
      RecvRequest[i] =  MPI::REQUEST_NULL;
    }else{    
      data_out[i].resize( cnt );
      RecvRequest[i] = MPI::COMM_WORLD.Irecv(&(data_out[i][0]), cnt, MPI::UNSIGNED, i, 13);
    }
  }
  
  // setup non-blocking sends
  vector<MPI::Request> SendRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_shared(i);

    if((i==MyRank)||(cnt==0)){
      SendRequest[i] =  MPI::REQUEST_NULL;
    }else{    
      SendRequest[i] = MPI::COMM_WORLD.Isend(&(data_in[i][0]), cnt, MPI::UNSIGNED, i, 13);
    }
  }
  
  MPI::Request::Waitall(NProcs, &(RecvRequest[0]));
  MPI::Request::Waitall(NProcs, &(SendRequest[0]));
  
}
// Updata all halo values stored in data_in
void Mesh::halo_update(const vector< vector<unsigned> >& data_in, const vector<unsigned>& data_cnt,
                       vector< vector<unsigned> >& data_out){
  
  data_out.clear();
  data_out.resize( NProcs );

  // Calculate the number of items being sent per node
  unsigned stride=0;
  for(unsigned i=0; i<(unsigned)NProcs; i++){
    if( data_in[i].empty() )
      continue;
    else{
      stride = data_in[i].size() / num_nodes_shared(i);
      break;
    }
  }
    
  // setup non-blocking receives
  vector<MPI::Request> RecvRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*data_cnt[i];
    
    if((i==MyRank)||(cnt==0)){
      RecvRequest[i] =  MPI::REQUEST_NULL;
    }else{
      data_out[i].resize(cnt);
      RecvRequest[i] = MPI::COMM_WORLD.Irecv(&(data_out[i][0]), cnt, MPI::UNSIGNED, i, 13);
    }

  }
  
  // setup non-blocking sends
  vector<MPI::Request> SendRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_shared(i);
    
    if((i==MyRank)||(cnt==0)){
      SendRequest[i] =  MPI::REQUEST_NULL;
    }else{
      SendRequest[i] = MPI::COMM_WORLD.Isend(&(data_in[i][0]), cnt, MPI::UNSIGNED, i, 13);
    }
  }
  
  MPI::Request::Waitall(NProcs, &(RecvRequest[0]));
  MPI::Request::Waitall(NProcs, &(SendRequest[0]));
  
}
void Mesh::halo_update(const vector< vector<samfloat_t> >& data_in, vector< vector<samfloat_t> >& data_out){
  
  data_out.clear();
  data_out.resize( NProcs );
  
  // Calculate the number of items being sent per node
  unsigned stride=0;
  for(unsigned i=0; i<(unsigned)NProcs; i++){
    if( data_in[i].empty() )
      continue;
    else{
      stride = data_in[i].size() / num_nodes_shared(i);
      break;
    }
  }
  
  // setup non-blocking receives
  vector<MPI::Request> RecvRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_halo(i);
    
    if((cnt==0)||(i==MyRank)){
      RecvRequest[i] =  MPI::REQUEST_NULL;
    }else{    
      data_out[i].resize( cnt );
      RecvRequest[i] = MPI::COMM_WORLD.Irecv(&(data_out[i][0]), cnt, MPI::SAMFLOAT, i, 13);
    }
  }
  
  // setup non-blocking sends
  vector<MPI::Request> SendRequest(NProcs);
  for(int i = 0; i<NProcs; i++){
    int cnt = stride*num_nodes_shared(i);
    
    if((cnt==0)||(i==MyRank)){
      SendRequest[i] =  MPI::REQUEST_NULL;
      continue;
    }else{
      SendRequest[i] = MPI::COMM_WORLD.Isend(&(data_in[i][0]), cnt, MPI::SAMFLOAT, i, 13);
    }
  }
  
  MPI::Request::Waitall(NProcs, &(RecvRequest[0]));
  MPI::Request::Waitall(NProcs, &(SendRequest[0]));
  
}

unsigned Mesh::min_node_owner(const unsigned elem){
  const vector<unn_t>& nodes = element_list[elem].get_enlist();
  unsigned len = nodes.size();
  assert(len>0);
  
  unsigned min = node_list[ unn2gnn(nodes[0]) ].get_current_owner();
  for(unsigned i=1; i<len; i++){
    unsigned T = node_list[ unn2gnn(nodes[i]) ].get_current_owner();
    min = min<T?min:T;
  }
  
  return min;
}

#define ALPHA_CONST 0.20412414523193154

#define mtetin_fc F77_FUNC(mtetin, MTETIN)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void mtetin_fc(samfloat_t *, samfloat_t *, samfloat_t*, samfloat_t*, 
	       samfloat_t*, samfloat_t*, samfloat_t*, samfloat_t*, samfloat_t*);

samfloat_t Mesh::element_functional(const unsigned eid){
  unsigned char flags = element_list[eid].get_flags();

  // element_functional only calculated for volume elements
  if( flags&ELM_SURFACE ){
    WARNING("element functional only calculated for volume elements");
    return 0.0;
  }

  ECHO("Calculating element function for element "<<eid);
  assert(eid<element_list.size());
  
  samfloat_t fxnl=0.0;
  if (dimension == 3)
  {
    // Use Adrians fortran subroutine
    samfloat_t X[4], Y[4], Z[4], VOL, AREAS[4], RADIUS, QUALTY, L[6];
    vector<samfloat_t> M(9, 0.0);
    
    vector<unn_t> nodes( element_list[eid].get_enlist() );
    assert( nodes.size() == 4 );
    for(unsigned i=0; i<4; i++)
      nodes[i] = unn2gnn( nodes[i] );
    
    // Get coordinates and average metric
    for(unsigned i=0; i<4; i++){
      X[i] =  node_list[ nodes[i] ].get_x();
      Y[i] =  node_list[ nodes[i] ].get_y();
      Z[i] =  node_list[ nodes[i] ].get_z();
      
      const vector<samfloat_t>& metric = node_list[ nodes[i] ].get_metric();
      
      CHECK( metric.size() );
      assert(metric.size() == 9);
      for(unsigned j=0; j<9; j++){
	M[j] += metric[j];
      }
    }
    
    for(unsigned j=0; j<9; j++)
      M[j] *= 0.25;
    
    // Adrians routine
    mtetin_fc(X, Y, Z, &(M[0]), &VOL, AREAS, L, &RADIUS, &QUALTY);
    CHECK(VOL);
    CHECK(RADIUS);
    CHECK(QUALTY);
    
    // From this the functional can be calculated as
    fxnl = 1.0;
    for(unsigned i=0; i<6; i++)
      fxnl += (1.0 - L[i])*(1.0 - L[i]);
    fxnl*=0.5;
    assert(RADIUS != 0);
    fxnl += (1.0 - ALPHA_CONST/RADIUS)*(1.0 - ALPHA_CONST/RADIUS);
  }
  else if (dimension == 2)
  {
    samfloat_t X[3], Y[3], M1[3], M2[3], M3[3];

    vector<unn_t> nodes( element_list[eid].get_enlist() );
    assert( nodes.size() == 3 );
    for(unsigned i=0; i<3; i++)
      nodes[i] = unn2gnn( nodes[i] );
    
    // Get coordinates and average metric
    for(unsigned i=0; i<3; i++)
    {
      X[i] =  node_list[ nodes[i] ].get_x();
      Y[i] =  node_list[ nodes[i] ].get_y();
    }
    const vector<samfloat_t>& metric1 = node_list[ nodes[0] ].get_metric();
    M1[0] = metric1[0]; M1[1] = metric1[1]; M1[2] = metric1[3];
    const vector<samfloat_t>& metric2 = node_list[ nodes[1] ].get_metric();
    M2[0] = metric2[0]; M2[1] = metric2[1]; M2[2] = metric2[3];
    const vector<samfloat_t>& metric3 = node_list[ nodes[2] ].get_metric();
    M3[0] = metric3[0]; M3[1] = metric3[1]; M3[2] = metric3[3];
    
    functional_2d f;
    fxnl = f.standard(X[0], Y[0], M1, X[1], Y[1], M2, X[2], Y[2], M3);
  }
  
  CHECK(fxnl);
  return fxnl;
}
#undef ALPHA_CONST

samfloat_t Mesh::elementVolume(const unsigned eid){
  // Warning...this only works for tets!!
  const std::vector<unn_t>& nodes = element_list[eid].get_enlist();
  samfloat_t vol;

  if (nodes.size() == 4)
  {
    const vector<samfloat_t>& r0 = node_list.unn( nodes[0] ).get_coord();
    const vector<samfloat_t>& r1 = node_list.unn( nodes[1] ).get_coord();
    const vector<samfloat_t>& r2 = node_list.unn( nodes[2] ).get_coord();
    const vector<samfloat_t>& r3 = node_list.unn( nodes[3] ).get_coord();
    
    samfloat_t a1 = r1[0]-r0[0];
    samfloat_t a2 = r1[1]-r0[1];
    samfloat_t a3 = r1[2]-r0[2];
    
    samfloat_t b1 = r2[0]-r0[0];
    samfloat_t b2 = r2[1]-r0[1];
    samfloat_t b3 = r2[2]-r0[2];
    
    samfloat_t c1 = r3[0]-r0[0];
    samfloat_t c2 = r3[1]-r0[1];
    samfloat_t c3 = r3[2]-r0[2];
    
    // volume = | r_a r_b r_c | / 6
    vol = (a1*(b2*c3 - b3*c2) - b1*(a2*c3 - a3*c2) + c1*(a2*b3 - a3*b2))/6.0;
  }
  else if (nodes.size() == 3)
  {
    const vector<samfloat_t>& r0 = node_list.unn( nodes[0] ).get_coord();
    const vector<samfloat_t>& r1 = node_list.unn( nodes[1] ).get_coord();
    const vector<samfloat_t>& r2 = node_list.unn( nodes[2] ).get_coord();
    samfloat_t a1 = r1[0]-r0[0];
    samfloat_t a2 = r1[1]-r0[1];
    
    samfloat_t b1 = r2[0]-r0[0];
    samfloat_t b2 = r2[1]-r0[1];

    vol = (a1*b2 - a2*b1) / 2.0;
  }
  else
  {
    ERROR("Only know how to compute the volume of triangles and tets. Sorry!");
  }
  
  // return the abs value as there may be orientation problems
  return fabs(vol);
}

samfloat_t Mesh::idealElementDensity(const unsigned eid)
{
  assert(element_list[eid].get_flags()&ELM_VOLUME);

  if (dimension == 3)
  {

    // First get the metric for this element. The metric for the element
    // is defined as average of the metric of all the nodes in that
    // element.
    samfloat_t a1 = 0.0;
    samfloat_t a2 = 0.0;
    samfloat_t a3 = 0.0;
    
    samfloat_t b1 = 0.0;
    samfloat_t b2 = 0.0;
    samfloat_t b3 = 0.0;
    
    samfloat_t c1 = 0.0;
    samfloat_t c2 = 0.0;
    samfloat_t c3 = 0.0;
    
    const std::vector<unn_t>& nodes = element_list[eid].get_enlist();  
    for(vector<unsigned>::const_iterator it=nodes.begin();it!=nodes.end();++it){
      const vector<samfloat_t>& metric = node_list.unn( *it ).get_metric();
      // it doesn't matter what order they are loaded in as since the
      // metric is symmetric
      a1+=metric[0];
      a2+=metric[1];
      a3+=metric[2];
      b1+=metric[3];
      b2+=metric[4];
      b3+=metric[5];
      c1+=metric[6];
      c2+=metric[7];
      c3+=metric[8];
    }
    unsigned len = nodes.size();
    a1/=len;
    a2/=len;
    a3/=len;
    b1/=len;
    b2/=len;
    b3/=len;
    c1/=len;
    c2/=len;
    c3/=len;
    
    // det = | r_a r_b r_c |
    samfloat_t det = a1*(b2*c3 - b3*c2) - b1*(a2*c3 - a3*c2) + c1*(a2*b3 - a3*b2);
    CHECK(det);
    
    return sqrt(det)*elementVolume(eid);
  }
  else if (dimension == 2)
  {
    samfloat_t a1 = 0.0, a2 = 0.0, b1 = 0.0, b2 = 0.0;
    const std::vector<unn_t>& nodes = element_list[eid].get_enlist();  
    for(vector<unsigned>::const_iterator it=nodes.begin();it!=nodes.end();++it)
    {
      const vector<samfloat_t>& metric = node_list.unn( *it ).get_metric();
      // it doesn't matter what order they are loaded in as since the
      // metric is symmetric
      a1+=metric[0];
      a2+=metric[1];
      b1+=metric[2];
      b2+=metric[3];
    }
    unsigned len = nodes.size();
    a1 /= len; a2 /= len; b1 /= len; b2 /= len;
    samfloat_t det = a1*b2 - a2*b1;
    return sqrt(det)*elementVolume(eid);
  }
}

// Returns a node->element list
deque< set<unsigned> > Mesh::mknelist(){
  deque< set<unsigned> > nelist( __num_nodes_total );
  
  // Loop through all the elements
  unsigned eid = 0;
  for(ElementVector<Element>::iterator ie = element_list.begin(); ie != element_list.end(); ++ie, ++eid){
    unsigned char flag = (*ie).get_flags();
    
    // Default behaviour only deals in volume elements
    if(flag&ELM_VOLUME){
      
      // Get node list
      const vector<unn_t>& nodes = (*ie).get_enlist();
      
      for(vector<unn_t>::const_iterator in = nodes.begin(); in != nodes.end(); ++in){
	unsigned gnn = unn2gnn(*in);
	assert(gnn<__num_nodes_total);
	nelist[ gnn ].insert( eid );
      }
    }
  }
  
  return nelist;
}

// Returns a node->element list in a map with an option to use unn's
map<unsigned, set<unsigned> > Mesh::mknelist(const bool use_unn){
  
  map<unsigned, set<unsigned> > nelist;
  
  // Loop through all the elements
  unsigned eid = 0;
  for(ElementVector<Element>::iterator ie = element_list.begin(); ie != element_list.end(); ++ie, ++eid){
    unsigned char flag = (*ie).get_flags();
    
    // Default behaviour only deals in volume elements
    if(flag&ELM_VOLUME){
      
      // Get node list
      const vector<unn_t>& nodes = (*ie).get_enlist();
      
      for(vector<unn_t>::const_iterator in = nodes.begin(); in != nodes.end(); ++in){
	if(use_unn)
	  nelist[ *in ].insert( eid );
	else
	  nelist[ unn2gnn(*in) ].insert( eid );
      }
    }
  }
  
  return( nelist );
}

//
// Returns a node->node list in a map an option to use unns
//
map<unsigned, set<unsigned> > Mesh::mknnlist(const bool use_unn){
  map<unsigned, set<unsigned> > edges;

  // Loop through elements
  for(ElementVector<Element>::iterator it = element_list.begin(); it != element_list.end(); ++it){
    vector<unn_t> enlist( (*it).get_enlist() );
    
    for(vector<unn_t>::iterator i = enlist.begin(); i != enlist.end(); ++i){
      unsigned gnn = unn2gnn(*i);
      
      if(node_list[gnn].get_current_owner() != MyRank)
	continue;
      
      for(vector<unn_t>::iterator j = enlist.begin(); j != enlist.end(); ++j){
	
	if( (*i) != (*j) ){
	  
	  if(use_unn){
	    edges[*i].insert( *j );
	  }else{
	    edges[gnn].insert( unn2gnn(*j) );
	  }
	}

      } // end of element. 

    }
  }

#ifndef NDEBUG
  ECHO("\t Look at node-node list");
  for(map<unsigned, set<unsigned> >::iterator it=edges.begin(); it!=edges.end(); ++it){
    ECHO((*it).first);
    for(set<unsigned>::iterator jt=(*it).second.begin(); jt!=(*it).second.end(); ++jt){
      ECHO(*jt);
    }
  }
#endif

  return edges;
}

//
// Function returns true if the mesh is consistant and false otherwise.
//
bool Mesh::mesh_consistant(){
#ifndef NDEBUG
  ECHO("Testing that the meshs are complete");
  
  unsigned elem = 0;
  for(ElementVector<Element>::const_iterator et = element_list.begin(); et!=element_list.end(); ++et){
    ECHO("Checking element "<<elem);
    
    // regular nodes
    ECHO("regular nodes: ");
    const vector<unn_t>& enl = (*et).get_enlist();
    for(vector<unn_t>::const_iterator it = enl.begin(); it != enl.end(); ++it){
      ECHO(*it);
      if(!node_list.contains_unn( *it ))
	return(false);
    }
    
    // Pressure nodes.
    ECHO("pressure nodes: ");
    const vector<unn_t>& enl2 = (*et).get_MFenlist();
    for(vector<unn_t>::const_iterator it = enl2.begin(); it != enl2.end(); ++it){
      ECHO(*it);      
      if(!MFnode_list.contains_unn( *it ))
	return(false);
    }
  }
  
  // regular halo
  ECHO("Testing regular halo");
  for(int p=0; p<NProcs; p++){
    for(set<unsigned>::const_iterator it = shared_nodes[p].begin(); it != shared_nodes[p].end(); ++it){
      if(!node_list.contains_unn(*it))
	return(false);
    }
    for(set<unsigned>::const_iterator it = halo_nodes[p].begin(); it != halo_nodes[p].end(); ++it){
      if(!node_list.contains_unn(*it))
	return(false);
    }
  }
  
  
  // pressure halo
  ECHO("Testing pressure halo");
  for(int p=0; p<NProcs; p++){
    for(set<unsigned>::const_iterator it = shared_pnodes[p].begin(); it != shared_pnodes[p].end(); ++it){
      if(!MFnode_list.contains_unn(*it))
	return(false);
    }
    for(set<unsigned>::const_iterator it = halo_pnodes[p].begin(); it != halo_pnodes[p].end(); ++it){
      if(!MFnode_list.contains_unn(*it))
	return(false);
    }
  }
#endif
  return(true);
}

void Mesh::set_functional_tol(samfloat_t tol){
  functional_tolerence = tol;
  return;
}

int Mesh::get_ncolga()
{
  int ncolga;

  ncolga = 0;
  for (int i = 0; i < NProcs; i++)
  {
    ncolga = ncolga + shared_nodes[i].size();
  }

  return ncolga;
}

int Mesh::get_nscate()
{
  int nscate;

  nscate = 0;

  for (int i = 0; i < NProcs; i++)
  {
    nscate = nscate + halo_nodes[i].size();
  }

  return nscate;
}


int Mesh::get_pncolga() const
{
  if(!mixed_formulation()){
    return 0;
  }

  int pncolga = 0;

  for(int i = 0;i < NProcs;i++)
  {
    pncolga = pncolga + shared_pnodes[i].size();
  }

  return pncolga;
}

int Mesh::get_pnscate() const
{
  if(!mixed_formulation()){
    return 0;
  }

  int pnscate = 0;

  for(int i = 0;i < NProcs;i++){
    pnscate = pnscate + halo_pnodes[i].size();
  }

  return pnscate;
}
