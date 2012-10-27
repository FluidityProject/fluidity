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
// provides interface for METIS
//

#include <vector>
#include <set>
#include <deque>
#include <assert.h>
#include <string.h>

#include "csr.h"
#include "samMetis.h"

using std::deque;
using std::vector;
using std::set;
using csr::Graph;

// METIS_PartGraphRecursive(int *, idxtype *, idxtype *, 
//                          idxtype *, idxtype *, int *, 
//                          int *, int *, int *, int *, idxtype *);
void Graph::PartGraphRecursive(int nparts, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Region growing
  options[3] = 1; // Early-Exit Boundary FM refinment
  options[4] = 0; // always set to 0

  METIS_PartGraphRecursive(&ncnt,  
			   bptr.begin(), edges.begin(), 
			   NULL, NULL, &wflag, &cflag,
			   &nparts, options,
			   &edgecut, noddom.begin());

}
 
// METIS_PartGraphKway(int *, idxtype *, idxtype *, 
//                     idxtype *, idxtype *, int *, 
//                     int *, int *, int *, int *, idxtype *);
void Graph::PartGraphKway(int nparts, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Multilevel recursive bisection
  options[3] = 3; /* Random boundary refinment that also minimizes 
		     the connectivity among the subdomains. */
  options[4] = 0; // always set to 0
  
  METIS_PartGraphKway(&ncnt,  
		      bptr.begin(), edges.begin(), 
		      NULL, NULL, &wflag, &cflag,
		      &nparts, options,
		      &edgecut, noddom.begin());

}

// METIS_PartGraphVKway(int *, idxtype *, idxtype *, 
//                      idxtype *, idxtype *, int *, 
//                      int *, int *, int *, int *, idxtype *);
void Graph::PartGraphVKway(int nparts, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Multilevel recursive bisection
  options[3] = 1; // Random boundary refinment
  options[4] = 0; // always set to 0
  
  METIS_PartGraphVKway(&ncnt,  
		       bptr.begin(), edges.begin(), 
		       NULL, NULL, &wflag, &cflag,
		       &nparts, options,
		       &edgecut, noddom.begin());

}

// METIS_mCPartGraphRecursive(int *, int *, idxtype *, idxtype *, 
//                            idxtype *, idxtype *, int *, 
//                            int *, int *, int *, int *, idxtype *);
void Graph::mCPartGraphRecursive(int nparts, int ncon, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 2; // Random
  options[3] = 1; // Early-Exit Boundary FM refinment
  options[4] = 0; // always set to 0

  METIS_mCPartGraphRecursive(&ncnt, &ncon,  
			     bptr.begin(), edges.begin(), 
			     NULL, NULL, &wflag, &cflag,
			     &nparts, options,
			     &edgecut, noddom.begin());

}

// METIS_mCPartGraphKway(int *, int *, idxtype *, idxtype *, 
//                       idxtype *, idxtype *, int *, 
//                       int *, int *, int *, int *, idxtype *);
void Graph::mCPartGraphKway(int nparts, int ncon, vector<float>& vcon, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Multilevel recursive bisection
  options[3] = 3; /* Random boundary refinment that also minimizes 
		     the connectivity among the subdomains. */
  options[4] = 0; // always set to 0
  
  METIS_mCPartGraphKway(&ncnt, &ncon,
			bptr.begin(), edges.begin(), 
			NULL, NULL, &wflag, &cflag,
			&nparts, vcon.begin(), options,
			&edgecut, noddom.begin());

}

// METIS_WPartGraphRecursive(int *, idxtype *, idxtype *, 
//                          idxtype *, idxtype *, int *, 
//                          int *, int *, float *,int *, int *, idxtype *);
void Graph::WPartGraphRecursive(int nparts, vector<float>& pweights, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Region growing
  options[3] = 1; // Early-Exit Boundary FM refinment
  options[4] = 0; // always set to 0

  METIS_WPartGraphRecursive(&ncnt,  
			    bptr.begin(), edges.begin(), 
			    NULL, NULL, &wflag, &cflag,
			    &nparts, pweights.begin(), options,
			    &edgecut, noddom.begin());

}

// METIS_WPartGraphKway(int *, idxtype *, idxtype *, 
//                     idxtype *, idxtype *, int *, 
//                     int *, int *, float *, int *, int *, idxtype *);
void Graph::WPartGraphKway(int nparts, const vector<float>& pweights, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Multilevel recursive bisection
  options[3] = 3; /* Random boundary refinment that also minimizes 
		     the connectivity among the subdomains. */
  options[4] = 0; // always set to 0
  
  METIS_WPartGraphKway(&ncnt,  
		       bptr.begin(), edges.begin(), 
		       NULL, NULL, &wflag, &cflag,
		       &nparts, (float *)(pweights.begin()), options,
		       &edgecut, noddom.begin());

}

// METIS_WPartGraphVKway(int *, idxtype *, idxtype *, 
//                       idxtype *, idxtype *, int *, 
//                       int *, int *, float *, int *, int *, idxtype *);
void Graph::WPartGraphVKway(int nparts, vector<float>& pweights, vector<int>& noddom){
  
  int wflag = 0, cflag = 0, edgecut;
  int options[5];
  
  int ncnt = nnodes;

  options[0] = 1; // Read remaining options
  options[1] = 3; // Sorted heavy edge matching
  options[2] = 1; // Multilevel recursive bisection
  options[3] = 1; // Random boundary refinment
  options[4] = 0; // always set to 0
  
  METIS_WPartGraphVKway(&ncnt,  
			bptr.begin(), edges.begin(), 
			NULL, NULL, &wflag, &cflag,
			&nparts, pweights.begin(), options,
			&edgecut, noddom.begin());

#ifndef METIS3
  METIS_WPartGraphVKway(&ncnt,  
			bptr.begin(), edges.begin(), 
			&nparts, pweights.begin(), options,
			&edgecut, noddom.begin());
#endif

}

// Sparse matrix reording Routines
void Graph::NodeND(vector<int>& perm, vector<int>& iperm){
  int cflag = 0;
  int options[8];
  
  int ncnt = nnodes;
  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  options[3] = 2;
  options[4] = 0;
  options[5] = 1;
  options[6] = 0;
  options[7] = 1;
  
  
  METIS_NodeND(&ncnt, bptr.begin(), edges.begin(), 
	       &cflag, options, perm.begin(), iperm.begin());
}














