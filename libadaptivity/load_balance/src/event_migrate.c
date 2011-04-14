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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include "migrate.h"
#include "mesh.h"
#include "LUTs.h"

/* This is the migrate routine to be called by EVENT. After the second 
 * route has been called, and EVENT has it's new mesh, it is responsible for
 * creating it's new node-element list.
 */
void emigrate1_(const int *PNodes, const int *NNodes, const int *NElems,  
		const int *scatter, const int *ATOREC, const int *NScatter,
		const int *gather, const int *ATOSEN, const int *NGather, 
		const int *ENLIST, const int *ENLBAS, const int *SZENLS,
		const int *NELIST, const int *NELBAS, const int *SZNELS,
		const int *ITYP, const int *IAVR, const int *IMAT, const int *ISOR,
		const int *ZeroOffset, const int *NODDOM,
		const coord_t *NODX, const coord_t *NODY, const coord_t *NODZ,
		const metric_t *METRIC, const field_t *FIELDS, const int *TOTUNS,
		/* return to fortran...*/
		int *F_NewPNodes, int *F_NewNNodes, int *F_NewSZENLS,
		int *F_NewNScatter, int *F_NewNGather,
		int *IERROR){
  
  int MyRank, NProcs;
  int i,j, *PNodesPerProc, *NNodesPerProc, err;
  long int *G_scatter, *G_gather;
  int *scatterCnt, *gatherCnt, *tmptr;
  long int refnum;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
  
  /* Convert data to horizons internal graph format.
   */
  if((err=event2internal((const int)*NNodes,
			 (const int)*PNodes,
			 (const int)*NElems,
			 (const int)*TOTUNS,
			 (const int *)ITYP, 
			 (const int *)IAVR, 
			 (const int *)IMAT, 
			 (const int *)ISOR,
			 (const coord_t *)NODX, 
			 (const coord_t *)NODY, 
			 (const coord_t *)NODZ,
			 (const metric_t *)METRIC, 
			 (const field_t *)FIELDS,
			 (const int *)ENLIST,
			 (const int *)ENLBAS,
			 (const int  )*SZENLS,
			 (const int *)NELIST,
			 (const int *)NELBAS,
			 (const int  )*SZNELS,
			 (const int *)ATOREC,
			 (const int *)scatter,
			 (const int)*NScatter,
			 (const int *)ATOSEN,
			 (const int *)gather,
			 (const int)*NGather,
			 (const int)*TOTUNS,
			 (const char)*ZeroOffset))<0){    
    *IERROR = err;
    return;
  }
  
  /* Send around new nodes according to NODDOM.
   */
  if( (err = migrate(NODDOM))<0){
    *IERROR = err;
    return;
  }
  
  /* Send back.
   */
  *F_NewPNodes   = NumNodes;
  *F_NewNNodes   = NumPNodes;
  *F_NewSZENLS   = NumElems;
  *F_NewNScatter = meshNSCAT;
  *F_NewNGather  = meshNGATH;

  return; 
}

void emigrate2_(const int *OldNNodes,
		int *scatter, int *ATOREC, int *NScatter,
		int *gather,  int *ATOSEN, int *NGather, 
		int *ENLIST,  int *ENLBAS, int *SZENLS,
		int *ITYP, int *IAVR, int *IMAT, int *ISOR,
		int *ZeroOffset, coord_t *NODX, coord_t *NODY, coord_t *NODZ,
		metric_t *METRIC, field_t *FIELDS,
		int *TOTUNS,
		int *F_NewPNodes, int *F_NewNNodes, int *F_NewSZENLS, 
		int *IERROR){
  int n,i,j,cnt, MyRank, NProcs;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);

  /* Write the new halo. */
  for(i=0;i<NProcs+1;i++) 
    ATOREC[i] = meshATOREC[i] + (*ZeroOffset);  
  for(i=0; i<meshNSCAT;i++) 
    scatter[i] = meshScatter[i];
  
  for(i=0;i<NProcs+1;i++)    
    ATOSEN[i] = meshATOSEN[i] + (*ZeroOffset);  
  for(i=0; i<meshNGATH;i++) 
    gather[i] = meshGather[i];
  
  /* Write node-wise values.
   */
  for(n=0; n<NumNodes; n++){
    NODX[n]            = node_list[n].x[0];
    NODY[n]            = node_list[n].x[1];
#ifndef MESH_2D
    NODZ[n]            = node_list[n].x[2];
#endif
    
#ifdef USING_ADAPTIVITY
#ifndef MESH_2D
    memcpy(METRIC+i*9, &node_list[i].metric[0], 9*sizeof(metric_t));
#else
    memcpy(METRIC+i*4, &node_list[i].metric[0], 4*sizeof(metric_t));
#endif
#endif
    
    memcpy(FIELDS+i*(*TOTUNS), node_list[n].fields, (*TOTUNS)*sizeof(field_t));
  }
  
  /* Wrap-up.
   */
  G2L_free();
  mesh_free();

  return; 
}
















