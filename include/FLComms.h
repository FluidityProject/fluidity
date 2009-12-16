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

#ifndef FLCOMMS_H
#define FLCOMMS_H

#include "confdefs.h"
#include "fmangle.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string.h>
#include <vector>

// Used for figure out the node mapping from tetra4 to tetra10
int fl_ilink2(int v0, int v1);

extern "C" {
  int ilink2_fc(const int *, const int*);
}

//               tag |            nodes per proc
typedef std::map<int, std::deque<std::vector<int> > > HaloDatabase_t;

#ifdef HAVE_MPI
//               tag |                   pattern        per proc |   datatype
typedef std::map<int, std::map<std::pair<int, int>, std::map<int, MPI::Datatype> > > MPIDatatypeCache_t;
#endif

//* Halo communication module
/** This class maintains a database of all halo's associated with
  different element discretisations on the base base mesh. It also
  provides functions for: creating additional halo types based on
  preregistered halos; updating solution nodes a halo index.
*/
class FLComms{
 public:
  //* Basic constructor
  FLComms();
  //* Destructer
  ~FLComms();

  //* Overloaded <<
  /** Prints out the internals of the class in a nice format.*/
  friend std::ostream &operator<<(std::ostream&, FLComms&);

  //* Clear the comms cache.
  /** Clear the comms cache. This should be used after mesh
      optimisation. */
  int ClearCache();

  //* Export a copy of the halo
  /** \param send vector of nodes to be sent
    * \param recv vector of nodes to be received
    * \return 0 is happy, otherwise very unhappy*/
  int ExportHalo(const int& tag, std::vector<std::vector<int> >& send, std::vector<std::vector<int> >& recv);

  //* Get the side of the halo
  /** This is generally used in conjunction with export_halo().
    \param tag halo identifier
    \param num_to_send total number of nodes to be sent
    \param num_to_recv total number of nodes to be received
    \return 0 is happy, otherwise very unhappy */
  int GetInfo(int tag, int *np, int *num_to_send, int *num_to_recv);

  int GetElementOwner(int tag, int eid);
  int GetNodeOwner(int tag, int nid);
  int GetNOwnedNodes(int tag);
  int GetGEN2UEN(int tag, const int *ENList, int NLocal, const int *gnn2unn, int *gen2uen);
  int GetGNN2UNN(int tag, int nnodes, int nnodp, int nblocks, int *gnn2unn);
  
  //* Get all current halo identifiers.
  /** \return All current halo identifers.
    */
  std::set<int> GetTags() const;
  //* Tests for existance of halo
  /** \param tag halo identifier
      \return true iff halo was been registered*/
  bool HaveTag(int tag) const;
  
  //* Initialise the module
  /** This member is necessary to allow an instance of the module to
    be initialized when declared as a global. This includes tasks that
    cannot be done by the constructer as it is called before
    MPI::Init()
    \return 0 is happy, otherwise very unhappy*/
  int Init();

  //* Returns true if the halo containers are empty.
  /** Returns true if the halo containers are empty. If IsParallel()
      then this returns true only if the partition under consideration
      doesn't have any mesh information.
   **/
  bool IsEmpty() const;

  //* Returns true if running in parallel
  /**  Returns true if running on two or more processors, otherwise
       false.
   **/
  bool IsParallel() const;

  int MergeSurfacePatches(int tag, int NSElements, int NLocal, const int *SENList, int *ids);

  int RegisterElements(int tag, int NNodes, int NElements, int NLocal, int *ENList);

  //* Register a new halo
  /** \param tag halo identifier
    * \param npnodes number of private nodes
    * \param send nodes to be sent
    * \param recv nodes to be received
    * \return 0 is happy, otherwise very unhappy*/
  int RegisterHalo(const int& tag, const unsigned int& npnodes,  const std::vector<std::vector<int> >& send, const std::vector<std::vector<int> >& recv);
  //* Register a new halo
  /** This assumes that the input is in the fortran friendly compressed format.
    \param tag halo identifier
    \param nprocs number of processors in computation
    \param send compressed list of node that to be sent
    \param bptr_send base pointers into send
    \param recv compressed list of node that to be received
    \param bptr_recv base pointers into recv
    \return 0 is happy, otherwise very unhappy*/
  int RegisterHalo(int tag, int NOwnedNodes,
		   const int *send, const int *bptr_send, 
		   const int *recv, const int *bptr_recv);
		   
  //* Unregister an existing halo.
  /** \param tag halo identifier
    * \return 0 is happy, otherwise very unhappy*/
  int UnregisterHalo(const int& tag);

  //* Reset halo information
  /** Delete all halos and cached information.
      \return 0 is happy, otherwise very unhappy*/
  int Reset();
  
  //* Reset halo information
  /** Delete halo and associated cached information.
    \param tag halo identifier
    \return 0 is happy, otherwise very unhappy*/
  int Reset(int tag);
  
  //* Get number of partitions.
  /** \return Number of partitions.
    */
  int GetNProcs() const;

  //* Set number of partitions.
  /** This sets the number of processes known by this halo. This
      should only be necessary to use when halo's are being
      manipulated in serial.
   */
  int SetNProcs(int NProcs);

  //* Test halo consistancy
  /** Stress tests halo for consistancy. This is useful when debugging
    a new halo. A possible weakness of this module is that it assumes
    that all nodes associated with the native partition are indexed
    before all others.
    \param tag halo identifier
    \param ref reference solution variable which is known to be consistant
    \param blk_size number of solution variables per node stored continiously
    \param NFields total number of fields
    \param stride stride between consective fields
    \param NNodes total number of nodes
    \param NPrivate number of nodes in native partition
    \return 0 is happy, otherwise very unhappy*/
  int Test(int tag, const void *ref, int blk_size, int NFields,
	   int stride, int NNodes, int NPrivate);
  
  int Test(int tag, const void *ref, int blk_size, int NFields,
	   int stride, int NNodes, int NPrivate, int NElements, const int *ENList);
  
  //* Generate a halo for Tetra10 given a halo for Tetra4
  /**
    \param tagT4 identifier for halo associated with Tetra4
    \param tetra4 element-node adjancy list for Tetra4 mesh
    \param tagT10 identifier for halo associated with Tetra10
    \param tetra10 element-node adjancy list for Tetra10 mesh
    \param NElements total number of elements in mesh
    \return 0 is happy, otherwise very unhappy*/
  int Tetra4ToTetra10(int tagT4, int NPrivateNodes, const int *tetra4, int tagT10, int *tetra10, int NElements);

  //* Update solution variables
  /** Use the halo to update the halo solution variables.
    \param tag halo identifier
    \param v array that comtains the solution variables
    \param blk_size number of solution variables per node stored continiously
    \param NFields total number of fields
    \param stride stride between consective fields
    \return 0 is happy, otherwise very unhappy*/
  int Update(int tag, void *v, int blk_size, int NFields, int stride);

  void VerboseOff();
  void VerboseOn();
  
 private:
  int GetNextCommTag();

  bool verbose, initialised;

  HaloDatabase_t send, recv;
  HaloDatabase_t element_send, element_recv;
#ifdef HAVE_MPI
  MPIDatatypeCache_t send_type, recv_type;
#endif
  
  int NProcs;
  size_t MyRank;
  int fl_comm_tag;

  std::map<int, int> NOwnedNodes;
  
  std::map<int, std::deque<int> > node_ownership;
  std::map<int, std::deque<int> > element_ownership;
  
   /** Generates and caches a MPI::Datatype to fit the pattern defined by the arguments:
    \param tag halo identifier
    \param blk_sizenumber of solution variables per node stored continiously
    \param NFields total number of fields
    \param stride stride between consective fields
    \return 0 is happy, otherwise very unhappy */
  int CreateTypes(int tag, int blk_size, int NFields, int stride);
};

//* Stores an edge using node id's.
/** This is a simple edge class useful when sorting edges. Directional
    information is not perserved this edge (a,b)==(b,a)
*/
class SimpleEdge: public std::pair<size_t, size_t>{
 public: 
  SimpleEdge(size_t, size_t);
  SimpleEdge(const SimpleEdge&);
  SimpleEdge(const std::pair<size_t, size_t>&);

  friend std::ostream &operator<<(std::ostream& out, const SimpleEdge& in);
  SimpleEdge &operator=(const SimpleEdge& in);
  SimpleEdge &operator=(const std::pair<size_t, size_t>& in);
  bool operator==(const SimpleEdge& in) const;
  bool operator!=(const SimpleEdge& in) const;
  bool operator<(const SimpleEdge& in) const;
};

//* Stores a four node tetrahedron using node id's.
/** This is a simple tetra 4 class useful when sorting elements.
*/
class SimpleTetra4{
 public: 
  SimpleTetra4(size_t, size_t, size_t, size_t);
  SimpleTetra4(const SimpleTetra4&);
  
  friend std::ostream &operator<<(std::ostream& out, const SimpleTetra4& in);
  SimpleTetra4 &operator=(const SimpleTetra4& in);
  bool operator==(const SimpleTetra4& in) const;
  bool operator!=(const SimpleTetra4& in) const;
  bool operator<(const SimpleTetra4& in) const;
  size_t nodes[4];
};

extern FLComms halo_manager;

#endif
