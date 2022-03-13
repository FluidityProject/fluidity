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
#include <confdefs.h>

#include "sam_mpi.h"
#include <cassert>
#include <cstdlib>
#include <vector>
#include <string>
#include <typeinfo>
#include <iostream>

#include "c++debug.h"

// In-place send-receive.
template<class T>
void allSendRecv(std::vector< std::vector<T> >& inout){
  int MyRank, NProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
  assert(NProcs == inout.size());
  
  // Inform all process how much data they are receiving from every
  // other process. This is implemented using MPI_Alltoall because
  // some implementations of MPI are a bit restrictive on the number
  // of non-blocking communications that can be initiated.
  std::vector<int> send_count(NProcs), recv_count(NProcs);
  for(size_t p=0;p<NProcs;p++){
    send_count[p] = inout[p].size();
  }
  MPI_Alltoall(&(send_count[0]), 1, MPI_INT,
               &(recv_count[0]), 1, MPI_INT, MPI_COMM_WORLD);

  // Receiving buffer space.
  std::vector< std::vector<T> > recvBuffer(NProcs);

  // MPI_Status sendStatus, recvStatus;
  std::vector<MPI_Request> sendRequest(NProcs, MPI_REQUEST_NULL);
  std::vector<MPI_Request> recvRequest(NProcs, MPI_REQUEST_NULL);
  MPI_Status recvStatus;

  // Determine the MPI datatype on the fly.
  MPI_Datatype mpi_type;
  if(typeid(T) == typeid(signed char)){
    mpi_type =  MPI_CHAR;
  }else if(typeid(T) == typeid(char)){
    mpi_type =  MPI_CHAR;
  }else if(typeid(T) == typeid(signed short int)){
    mpi_type =  MPI_SHORT;
  }else if(typeid(T) == typeid(signed int)){
    mpi_type =  MPI_INT;
  }else if(typeid(T) == typeid(int)){
    mpi_type =  MPI_INT;
  }else if(typeid(T) == typeid(signed long int)){
    mpi_type =  MPI_LONG;
  }else if(typeid(T) == typeid(long)){
    mpi_type =  MPI_LONG;
  }else if(typeid(T) == typeid(unsigned char)){
    mpi_type =  MPI_UNSIGNED_CHAR;
  }else if(typeid(T) == typeid(unsigned short int)){
    mpi_type =  MPI_UNSIGNED_SHORT;
  }else if(typeid(T) == typeid(unsigned int)){
    mpi_type =  MPI_UNSIGNED;
  }else if(typeid(T) == typeid(unsigned)){
    mpi_type =  MPI_UNSIGNED;
  }else if(typeid(T) == typeid(unsigned long int)){
    mpi_type =  MPI_UNSIGNED_LONG;
  }else if(typeid(T) == typeid(float)){
    mpi_type =  MPI_FLOAT;
  }else if(typeid(T) == typeid(double)){
    mpi_type =  MPI_DOUBLE;
  }else if(typeid(T) == typeid(long double)){
    mpi_type =  MPI_LONG_DOUBLE;
  }else{
    std::cerr << "ERROR: illegal type " << typeid(T).name()
              << " passed into allSendRecv(std::vector< std::vector<T> >& inout)"
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  const MPI_Datatype const_mpi_type = mpi_type;
  
  // Send
  for(size_t p=0;p<NProcs;p++){
    if(send_count[p]){
      MPI_Isend(&(inout[p][0]), send_count[p], const_mpi_type, p, 1, MPI_COMM_WORLD, &(sendRequest[p]));
    }
  }
  
  // Receive
  for(size_t p=0;p<NProcs;p++){
    if(recv_count[p]){
      // Allocate receive space.
      recvBuffer[p].resize(recv_count[p]);
      
      // initiate non-blocking receive
      MPI_Irecv(&(recvBuffer[p][0]), recv_count[p], const_mpi_type, p, 1, MPI_COMM_WORLD, &(recvRequest[p]));
    }
  }
  MPI_Waitall(NProcs, &(sendRequest[0]), MPI_STATUSES_IGNORE); // Wait all sends
  MPI_Waitall(NProcs, &(recvRequest[0]), MPI_STATUSES_IGNORE); // Wait all receives

  inout.swap(recvBuffer);
 
  return;
}

// In-place send-receive.
template<class T>
void allSendRecv(std::vector< std::set<T> >& inout){
  unsigned len = inout.size();

  // Write to vectors
  std::vector< std::vector<T> > in(len);
  for(unsigned i=0; i<len; i++){
    in[i].resize( inout[i].size() );
    unsigned j=0;
    for(typename std::set<T>::const_iterator it=inout[i].begin(); it!=inout[i].end(); ++it){
      in[i][j++] = *it;
    }
    inout[i].clear();
  }
  
  // Communicate
  allSendRecv(in);
  
  // Write to sets
  for(unsigned i=0; i<len; i++){
    for(typename std::vector<T>::const_iterator it=in[i].begin(); it!=in[i].end(); ++it){
      inout[i].insert(*it);
    }
    in[i].clear();
  }
  
  return;
}

template<class T>
void allSendRecv(const std::vector< std::vector<T> >& sendBuffer, 
		 std::vector< std::vector<T> >& recvBuffer){
  int MyRank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  size_t NProcs = sendBuffer.size();
  recvBuffer.clear();
  recvBuffer.resize(NProcs);
  
  // Determine the MPI datatype on the fly.
  MPI_Datatype mpi_type;

  if(typeid(T) == typeid(signed char)){
    mpi_type =  MPI_CHAR;
  }else if(typeid(T) == typeid(char)){
    mpi_type =  MPI_CHAR;
  }else if(typeid(T) == typeid(signed short int)){
    mpi_type =  MPI_SHORT;
  }else if(typeid(T) == typeid(signed int)){
    mpi_type =  MPI_INT;
  }else if(typeid(T) == typeid(int)){
    mpi_type =  MPI_INT;
  }else if(typeid(T) == typeid(signed long int)){
    mpi_type =  MPI_LONG;
  }else if(typeid(T) == typeid(long)){
    mpi_type =  MPI_LONG;
  }else if(typeid(T) == typeid(unsigned char)){
    mpi_type =  MPI_UNSIGNED_CHAR;
  }else if(typeid(T) == typeid(unsigned short int)){
    mpi_type =  MPI_UNSIGNED_SHORT;
  }else if(typeid(T) == typeid(unsigned int)){
    mpi_type =  MPI_UNSIGNED;
  }else if(typeid(T) == typeid(unsigned)){
    mpi_type =  MPI_UNSIGNED;
  }else if(typeid(T) == typeid(unsigned long int)){
    mpi_type =  MPI_UNSIGNED_LONG;
  }else if(typeid(T) == typeid(float)){
    mpi_type =  MPI_FLOAT;
  }else if(typeid(T) == typeid(double)){
    mpi_type =  MPI_DOUBLE;
  }else if(typeid(T) == typeid(long double)){
    mpi_type =  MPI_LONG_DOUBLE;
  }else{
    std::cerr << "ERROR: illegal type " << typeid(T).name()
              << " passed into allSendRecv(std::vector< std::vector<T> >& inout)"
              << std::endl;
    exit(-1);
  }  
  const MPI_Datatype const_mpi_type = mpi_type;

  // Inform all process how much data they are receiving from every
  // other process. This is implemented using MPI_Alltoall because
  // some implementations of MPI are a bit restrictive on the number
  // of non-blocking communications that can be initiated.
  std::vector<int> send_count(NProcs), recv_count(NProcs);
  for(size_t p=0;p<NProcs;p++){
    send_count[p] = sendBuffer[p].size();
  }
  MPI_Alltoall(&(send_count[0]), 1, MPI_INT,
                           &(recv_count[0]), 1, MPI_INT, MPI_COMM_WORLD);
  
  // MPI_Status sendStatus, recvStatus;
  std::vector<MPI_Request> sendRequest(NProcs, MPI_REQUEST_NULL);
  std::vector<MPI_Request> recvRequest(NProcs, MPI_REQUEST_NULL);
  MPI_Status recvStatus;
  
  // Send
  for(size_t p=0;p<NProcs;p++){
    if(send_count[p]){
      MPI_Isend(&(sendBuffer[p][0]), send_count[p], const_mpi_type, p, 1,MPI_COMM_WORLD, &(sendRequest[p]));
    }
  }

  // Receive
  for(size_t p=0;p<NProcs;p++){
    if(recv_count[p]){
      // Allocate receive space.
      recvBuffer[p].resize(recv_count[p]);
      
      // initiate non-blocking receive
      MPI_Irecv(&(recvBuffer[p][0]), recv_count[p], const_mpi_type, p, 1, MPI_COMM_WORLD, &(recvRequest[p]));
    }
  }

  MPI_Waitall(NProcs, &(sendRequest[0]), MPI_STATUSES_IGNORE); // Wait all sends
  MPI_Waitall(NProcs, &(recvRequest[0]), MPI_STATUSES_IGNORE); // Wait all receives
  return;
}

template<class T>
void allSendRecv(const std::vector< std::set<T> >& in, std::vector< std::set<T> >& out){
  unsigned len = in.size();
  out.clear();
  
  // Write to vectors
  std::vector< std::vector<T> > inout(len);
  for(unsigned i=0; i<len; i++){
    inout[i].resize( in[i].size() );
    unsigned j=0;
    for(typename std::set<T>::const_iterator it=in[i].begin(); it!=in[i].end(); ++it){
      inout[i][j++] = *it;
    }
  }
  
  // Communicate
  allSendRecv(inout);
  
  // Write to sets
  out.resize(len);
  for(unsigned i=0; i<len; i++){
    for(typename std::vector<T>::const_iterator it=inout[i].begin(); it!=inout[i].end(); ++it){
      out[i].insert(*it);
    }
    inout[i].clear();
  }
  
  return;
}

#undef TAG_33
#undef TAG_44
