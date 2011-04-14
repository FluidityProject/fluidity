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
/* Function to convert the mesh from FLUIDITY's matrix format
   to horizons. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<mpi.h>
#include"mesh.h"

#define TAG_HALO 1

int event2internal(const int NNOD, const int PNOD, const int NELM,
		   const int field_len, 
		   const int *ITYP, const int *IAVR, const int *IMAT, const int *ISOR,
		   const coord_t *X, const coord_t *Y, const coord_t *Z, 
		   const metric_t *metric, const field_t *field,
		   const int *ENLIST, const int *ENLBAS, const int SZENLS,
		   const int *NELIST, const int *NELBAS, const int SZNELS,
		   const int *ATOREC, const int *SCATER, const int NSCAT,
		   const int *ATOSEN, const int *GATHER, const int NGATH,
		   const int TOTUNS, const char ZeroOffset){
  
  int i, j, cnt;
  int MyRank, NProcs, noffset, *NperP;
  MPI_Request *hSendRequest=0, *hRecvRequest=0;
  MPI_Status  *hStatus=0;

  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);

  if( (hRecvRequest = (MPI_Request *)malloc(NProcs*sizeof(MPI_Request)))==NULL){
    perror("malloc");
    return(-1);
  }
  if( (hStatus = (MPI_Status *)malloc(NProcs*sizeof(MPI_Status)))==NULL){
    perror("malloc");
    return(-12);
  }

  /***********************************************************
   * Set up global node
   */
  if( (NperP = (int *)malloc(NProcs*sizeof(int))) == NULL ){
    perror("malloc");
    return(-1);
  }
  MPI_Allgather((int *)&PNOD, 1, MPI_INT, NperP, 1, MPI_INT, MPI_COMM_WORLD);

  noffset=0;
  for(i=0;i<MyRank;i++)
    noffset+=NperP[i];

  /******************************/
  /* convert node information ***/
  /******************************/
  NumNodes  = NNOD;
  NumPNodes = PNOD;
  NumElems  = NELM;
  
  /* Allocate space for the node list.
   */
  if((node_list = (node_t *)malloc(NNOD*sizeof(node_t)))==NULL){
    perror("fluidity2internal: malloc");
    return(-1);
  }
  node_list_len = NNOD;

  /* Write node data. 
   */
  for(i=0;i<NNOD;i++){

    /* global node numbering and ownership*/
    if(i<PNOD){
      node_list[i].UNN = noffset + i;
      node_list[i].owner = MyRank;
    }

    node_list[i].NCE = NELBAS[i+1]-NELBAS[i];
    if( (node_list[i].elems = 
	 (int *)malloc(node_list[i].NCE*sizeof(int)))==NULL ){
      perror("malloc");
      return(-2);
    }
    memcpy(node_list[i].elems, NELIST+NELBAS[i]-ZeroOffset, 
	   node_list[i].NCE*sizeof(int));
        
    /* Nodewise field values.
     */
    node_list[i].nfields = TOTUNS;
    if( (node_list[i].fields=
	 (field_t *)malloc(node_list[i].nfields*sizeof(field_t)))==NULL ){
      perror("malloc");
      return(-3);
    }
    memcpy(node_list[i].fields, field+i*field_len, 
	   field_len*sizeof(field_t));
    
    node_list[i].x[0] = X[i];
    node_list[i].x[1] = Y[i];
#ifndef MESH_2D
    node_list[i].x[2] = Z[i];
#endif
    
    /* Nodewise metric used for adaptivity
     */
#ifdef USING_ADAPTIVITY
    memcpy(&node_list[i].METRIC[0], metric+i*9, 
	   9*sizeof(metric_t));
#endif
  }

  /* Write owner information.
   */
  for(i=0;i<NProcs;i++)
    for(j=ATOREC[i]-ZeroOffset;j<ATOREC[i+1]-ZeroOffset;j++)
      node_list[SCATER[j]-ZeroOffset].owner=i;
    
  /* Write halo. */
  meshNGATH = NGATH;
  if( (meshATOSEN = (int *)realloc(meshATOSEN, (NProcs+1)*sizeof(int)))==NULL){
    perror("realloc");
    return(-4);
  }
  if( (meshGather= (unn_t *)realloc(meshGather, meshNGATH*sizeof(unn_t)))==NULL ){
    perror("realloc");
    return(-5);
  }
  for(i=0;i<=NProcs;i++)
    meshATOSEN[i] = ATOSEN[i]-ZeroOffset;
  for(i=0;i<meshNGATH;i++)
    meshGather[i] = GATHER[i]+noffset;

  meshNSCAT = NSCAT;
  if( (meshATOSEN = (int *)realloc(meshATOSEN, (NProcs+1)*sizeof(int)))==NULL){
    perror("realloc");
    return(-6);
  }
  if( (meshScatter= (unn_t *)realloc(meshScatter, meshNGATH*sizeof(unn_t)))==NULL ){
    perror("realloc");
    return(-7);
  }
  for(i=0;i<=NProcs;i++)
    meshATOSEN[i] = ATOREC[i]-ZeroOffset;
  for(i=0;i<NProcs;i++){                /* GNN of scatter */
    cnt = meshATOSEN[i+1]-meshATOSEN[i];
    if(cnt>0){
      MPI_Irecv(meshGather+meshATOSEN[i], cnt, MPI_LONG, i, TAG_HALO, 
		MPI_COMM_WORLD, hRecvRequest+i);
    }else{
      hRecvRequest[i] = MPI_REQUEST_NULL;
    }
  }
  for(i=0;i<NProcs;i++){
    cnt = meshATOSEN[i+1]-meshATOSEN[i];
    if(cnt>0){
      MPI_Isend(meshScatter+meshATOSEN[i], cnt, MPI_LONG, i, TAG_HALO, 
		MPI_COMM_WORLD, hSendRequest+i);
    }else{
      hSendRequest[i] = MPI_REQUEST_NULL;
    }
  }
  MPI_Waitall(NProcs, hSendRequest, hStatus);
  MPI_Waitall(NProcs, hRecvRequest,  hStatus);

  /*********************************/
  /* convert element information ***/
  /*********************************/
  if((element_list = (element_t *)malloc(NELM*sizeof(element_t))) == NULL){
    perror("fluidity2internal: malloc");
    return(-4);
  }
  element_list_len = NELM;
  
  for(i=0;i<NELM;i++){
        
    /* Number of nodes in element */
    element_list[i].NodeCnt = 
      ENLBAS[i+1]-ENLBAS[i];

    element_list[i].property[0] = ITYP[i];
    element_list[i].property[1] = IAVR[i]; /* Element region (material) 
				       information. AKA ELMREG */
    element_list[i].property[2] = IMAT[i];
    element_list[i].property[3] = ISOR[i];
    
#ifdef USING_DISCONTINUOUS
    /* We need to fill in these values in teh case of discontinuous methods.
       int nfields;
       field_t *fields;
    */
#endif
    /* A pointer into ElemNodeLIST indicating 
       the start of the node list for this element. */
    if( (element_list[i].nodes = 
	 (unn_t *)malloc(element_list[i].NodeCnt*sizeof(unn_t)))==NULL){
      perror("malloc");
      return(-3);
    }
    memcpy(element_list[i].nodes, ENLIST+ENLBAS[i]-ZeroOffset, 
	   element_list[i].NodeCnt*sizeof(int));
  }

  /* All appears to be well.
   */
  free(hSendRequest);
  free(hRecvRequest);
  free(hStatus);

  return(0);
}
