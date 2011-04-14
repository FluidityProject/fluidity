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
#ifndef H_MESH
#define H_MESH

/* This the the mesh stored in an internal form 
   designed to be more benificial for graph manipulation. 
   
   As can be seen below, the mesh is built using two classes,
   Node and Element. */
#include "confdefs.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <iostream>

#include "NodeVector.h"
#include "ElementVector.h"
#include "Graph.h"
#include "samtypes.h"
#include "Node.h"
#include "PressureNode.h"
#include "Element.h"

class Mesh{
 public:
  // Constructer & distructer
  Mesh();
  ~Mesh();
  
  // Overloaded operators

  //
  // Member functions.
  //
  
  // With no argument it returns true if we're using mixed
  // formulation. Otherwise it can be used to set the state.
  bool mixed_formulation() const;
  bool mixed_formulation(const bool);
  
  unsigned               max_nodepack_size() const;       // Estimate the size the the biggest node.
  unsigned               max_elementpack_size() const;    // Estimate the size the the biggest element.
  unsigned               max_pressurepack_size() const;   // Estimate the size the the biggest pressure.
  
  void                   add_node(const Node& );                        // Add a node to the end of the list.
  void                   add_node(const PressureNode& );
  void                   add_element(Element&);                         // Add an element to the end of the list.
  
  const Node&            get_node(const unsigned n) const;                 // Return a copy of the specified node
  const PressureNode&    get_pnode(const unsigned n) const;                // Return a copy of the specified pnode
  unsigned               num_nodes(const std::string);                     // {native, alien, total}
  std::vector<unn_t>     get_enlist(const unsigned);                       // returns a vector of nodes in some element
  unsigned               get_node_current_owner(const unsigned n);
  unsigned               get_node_future_owner(const unsigned n);
  unn_t                  get_node_unn(const unsigned);

  Element&               get_element(const unsigned e);
  unsigned               num_elements(const std::string& );             // {surface, volume, total}

  void                   add2shared_nodes(const unsigned short, const unn_t);
  void                   add2halo_nodes(const unsigned short, const unn_t);

  gnn_t                  unn2gnn(const unn_t unn);    // Given unn, what's gnn
  gnn_t                  MFunn2gnn(const unn_t unn);  // Given a MF-unn, what's gnn: Note..a const can't be appended
  
  // returned the minimum node owner of some  element
  unsigned min_node_owner(const unsigned elem);


  unsigned               num_nodes_shared();                        // Number of nodes shared with other processors
  unsigned               num_nodes_shared(const unsigned);          // Number of nodes shared with some processor.
  unsigned               num_nodes_halo();                          // Number of nodes in halo (from other processors)
  unsigned               num_nodes_halo(const unsigned);            // Number of nodes in halo (from some processor)
  unsigned               num_nodes_in_rank(const unsigned);         // Number of nodes on some rank
  void                   find_nodes_per_rank(const unsigned lcnt);  // Set the nodes_per_rank array taking lcnt as the local cnt

  std::set<int>   current_element_domain_set(const unsigned);// What domains currently use element x
  std::set<int>   future_element_domain_set(const unsigned); // What domains will use element x

  std::map<unn_t, int>   current_enode2domain_map(const unsigned);  // current mapping of the nodes in element x to domains
  std::map<unn_t, int>   future_enode2domain_map(const unsigned);   // future mapping of the nodes in element x to domains
  
  void  invent_pressure_mesh();
  void  formHalo2();

  samfloat_t element_functional(const unsigned);
  samfloat_t elementVolume(const unsigned);
  samfloat_t idealElementDensity(const unsigned);

  // Serious methods!
  void import_fluidity(const int, const int, const int, const int, const int,
		       const int,
		       const samfloat_t [],  const samfloat_t [], const samfloat_t [],
		       const samfloat_t [], const samfloat_t [],
		       const int [], const int,
		       const int [], const int [], const int,
		       const int [], const int [],
		       const int [], const int []);
  void import_pressure(const int, const int,
		       const samfloat_t [], const int,
		       const int [], const int [],
		       const int [], const int [],
		       const int [], const int []);
  void export_fluidity(int [],        int, int&, 
		       samfloat_t [], int, int&,
		       const int,
		       int&, int&, int&, int&, 
		       int&, int [],
		       int&, int [],
		       int&, int&, int&,
		       // Pressure stuff
		       int&, int&,
		       int&, int [],
		       int&, int [],
		       int&,
		       int&, int&, int&, samfloat_t []);
		       
  void export_halo(int* colgat, int* atosen, int* scater, int* atorec, const int* ncolga, const int* nscate, const int* nprocs);
  void export_phalo(int* pcolgat, int* patosen, int* pscater, int* patorec, const int* pncolga, const int* pnscate, const int* nprocs);
  
  void halo_update(const std::vector< std::vector<int> >&,       std::vector< std::vector<int> >&      );
  void halo_update(const std::vector< std::vector<unsigned> >&,  std::vector< std::vector<unsigned> >& );
  void halo_update(const std::vector< std::vector<unsigned> >&,  const std::vector<unsigned>&, std::vector< std::vector<unsigned> >& );
  
  void halo_update(const std::vector< std::vector<samfloat_t> >&,     
		   std::vector< std::vector<samfloat_t> >&    );

  std::vector<int> decomp(const std::vector<int>&);  
  void migrate(const std::vector<int> noddom);

  void set_functional_tol(samfloat_t);

  // Storage for mesh.
  NodeVector<Node>         node_list;
  ElementVector<Element>   element_list;
  NodeVector<PressureNode> MFnode_list;

  // fun methods
  int writeVTK(const std::string);

  void do_element_headcount();
  int get_ncolga();
  int get_nscate();
  int get_pncolga() const;
  int get_pnscate() const;
 private:
  // Mixed formulation mesh
  bool __mixed_formulation;
  bool mesh_consistant();
  samfloat_t functional_tolerence;

  // Node class tool-kit
  void                 set_node_unn(const unsigned, const unn_t);
  void                 set_node_current_owner(const unsigned n, const unsigned short o);
  void                 set_node_future_owner(const unsigned n, const unsigned short o);
  std::vector<samfloat_t> get_node_fields(const unsigned n);
  unsigned             get_enlist_size(const unsigned);
  void                 set_node_flags(const unsigned n, const unsigned char f);
  unsigned char        get_node_flags(const unsigned n);
  
  std::deque< std::set<unsigned> > mknelist();
  std::map<unsigned, std::set<unsigned> > mknelist(const bool use_unn);
  std::map<unsigned, std::set<unsigned> > mknnlist(const bool use_unn);
  
  int MyRank;
  int NProcs;

  // Data that defines a mesh. The reason a deque is choosen is because
  // there can be important resizes carried out. The meshes can get large
  // thus a realloc might wreak havoc on the computer.
//  std::deque<Node> node_list;  
  void do_node_headcount();
  unsigned __num_nodes_native;
  unsigned __num_nodes_alien;
  unsigned __num_nodes_total;

//  std::deque<Element> element_list;
  unsigned __num_elements_surface;
  unsigned __num_elements_volume;
  unsigned __num_elements_total;
  
  std::vector< std::set<unsigned> > shared_nodes;
  std::vector< std::set<unsigned> > halo_nodes;
  
  std::vector< std::set<unsigned> > shared_pnodes;
  std::vector< std::set<unsigned> > halo_pnodes;

  std::vector<int> nodes_per_rank;
  std::vector<int> unn_offsets;

  // Private functions.
  void calculate_submeshes( std::vector<Mesh>& );
  void fixate_elements(const unsigned, std::map<unsigned, Node>&, std::vector< std::map<unsigned, Node> >&);
  void fixate_nodes(std::vector< std::map<unsigned, Node> >&);
  void fixate_pressure();

  int dimension;
};

#endif
















