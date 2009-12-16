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

#include "FLComms.h"

using namespace std;

int fl_ilink2(int v0, int v1){
  return ilink2_fc(&v0, &v1);
}

#undef NDEBUG

void assert_no_pending_communication(const char *file, int line){
#ifdef DEBUG_FLCOMMS
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    MPI::Status iprobe_status;
    if(MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, iprobe_status)){
      cerr<<"ERROR - "<<file<<", "<<line<<": Found a pending communication - "
          <<"source="<<iprobe_status.Get_source()
          <<", int count="<<iprobe_status.Get_count(MPI::INT)
          <<", tag="<<iprobe_status.Get_tag()<<endl;
      MPI::COMM_WORLD.Abort(MPI::ERR_OTHER);
    }
    MPI::COMM_WORLD.Barrier();
  }
#endif
#endif
}

FLComms::FLComms(){
  NProcs = -1;
  MyRank = 0;
  verbose = false;
  initialised = false;
  fl_comm_tag = 0;
}

FLComms::~FLComms(){
}

int FLComms::ClearCache(){
#ifdef HAVE_MPI
  send_type.clear();
  recv_type.clear();
#endif
  return 0;
}

int FLComms::ExportHalo(const int& tag, vector<vector<int> >& send, vector<vector<int> >& recv){
  if(verbose)
    cout<<"int FLComms::ExportHalo(const int& tag, vector<vector<int> >& send, vector<vector<int> >& recv)\n";
    
  assert_no_pending_communication(__FILE__, __LINE__);
  
  if(!HaveTag(tag)){
    return -1;
  }

  send.clear();  send.resize(this->send[tag].size());
  for(size_t i=0;i<this->send[tag].size();i++){
    send[i]=this->send[tag][i];
  }
  recv.clear();  recv.resize(this->recv[tag].size());
  for(size_t i=0;i<this->recv[tag].size();i++){
    recv[i]=this->recv[tag][i];
  }
  
  return 0;
}

int FLComms::GetInfo(int tag, int *np, int *num_to_send, int *num_to_recv){
  if(verbose)
    cout<<"int FLComms::GetInfo(int tag, int *num_to_send, int *num_to_recv) const\n";

  Init();

  assert_no_pending_communication(__FILE__, __LINE__);

  *np = NProcs;  
  *num_to_send = 0;
  if(HaveTag(tag))
    for(size_t p=0;p<(size_t)NProcs;p++)
      *num_to_send += send[tag][p].size();
  
  *num_to_recv = 0;
  if(HaveTag(tag))
    for(size_t p=0;p<(size_t)NProcs;p++)
      *num_to_recv += recv[tag][p].size();

  return 0;
}

set<int> FLComms::GetTags() const{
  set<int> tags;
  for(HaloDatabase_t::const_iterator iter = send.begin();iter != send.end();iter++){
    if(HaveTag(iter->first)){
      tags.insert(iter->first);
    }
  }

  return tags;
}

bool FLComms::HaveTag(int tag) const{
  if(send.count(tag))
    if(recv.count(tag))
      return true;
  return false;
}

int FLComms::Init(){
  if(initialised)
    return 0;
  
  if(verbose)
    cout<<"int FLComms::Init()\n";
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    if(NProcs==-1){
      // The only reason this should be not -1 is if
      // it's value was set for serial
      // decomposition. In which case it shouldn't
      // be here.
      NProcs = MPI::COMM_WORLD.Get_size();
    }
    MyRank = MPI::COMM_WORLD.Get_rank();
  }else{
    if(NProcs==-1)
      NProcs = 1;
    MyRank = 0;
  }
#else
  NProcs = 1;
  MyRank = 0;
#endif
  
  initialised = true;
  return 0;
}

bool FLComms::IsEmpty() const{
  if(verbose)
    cout<<"bool FLComms::IsEmpty() const\n";
  
  for(map<int, int>::const_iterator it=NOwnedNodes.begin();it!=NOwnedNodes.end();++it)
    if(it->second>0)
      return false;
  
  return true;
}

bool FLComms::IsParallel() const{
  if(verbose)
    cout<<"bool FLComms::IsParallel() const\n";
  
#ifdef HAVE_MPI
  if(MPI::Is_initialized()){
    return (MPI::COMM_WORLD.Get_size()>1);
  }else{
    return false;
  }
#else
  return false;
#endif
}

void print_backtrace();

int FLComms::CreateTypes(int tag, int blk_size, int NFields, int stride){
  if(verbose)
    cout<<"int FLComms::CreateTypes(int tag="<<tag<<", int blk_size="<<blk_size
        <<", int NFields="<<NFields<<", int stride="<<stride<<")\n";

  assert(HaveTag(tag));
  assert_no_pending_communication(__FILE__, __LINE__);

#ifdef HAVE_MPI
  pair<int, int> pattern(blk_size, NFields);
  
  // Return right away if types for this pattern are already created
  if(send_type[tag].count(pattern)){
    if(verbose)
      cout<<"Have type\n";
    return 0;
  }

  if(verbose)
    cout<<"Creating new type: tag="<<tag<<", blk_size="<<blk_size
        <<", NFields="<<NFields<<", stride="<<stride<<endl;

  // Send types
  for(size_t p=0;p<(size_t)NProcs;p++){
    int hlen = send[tag][p].size();
    if(hlen==0)
      continue;
    
    // The length of each block being sent per node. This is constant.    
    vector<int> blens;
    
    // Offsets (displacements) to each block
    vector<int> disp;
    
    for(int i=0;i<NFields;i++){
      for(int j=0;j<hlen;j++){
        blens.push_back(blk_size);
        disp.push_back(i*stride + (send[tag][p][j]-1)*blk_size);
      }
    }

#ifdef DOUBLEP
    send_type[tag][pattern][p] = MPI::DOUBLE.Create_indexed(blens.size(), &(blens[0]), &(disp[0]));
#else
    send_type[tag][pattern][p] = MPI::FLOAT.Create_indexed(blens.size(), &(blens[0]), &(disp[0]));
#endif
    send_type[tag][pattern][p].Commit();
  }
  
  // Recv types
  for(size_t p=0;p<(size_t)NProcs;p++){
    int hlen = recv[tag][p].size();
    if(hlen==0)
      continue;
    
    // The length of each block being sent per node. This is constant.    
    vector<int> blens;
    
    // Offsets (displacements) to each block
    vector<int> disp;
    
    for(int i=0;i<NFields;i++){
      for(int j=0;j<hlen;j++){
        blens.push_back(blk_size);
        disp.push_back(i*stride + (recv[tag][p][j]-1)*blk_size);
      }
    }
    
#ifdef DOUBLEP
    recv_type[tag][pattern][p] = MPI::DOUBLE.Create_indexed(blens.size(), &(blens[0]), &(disp[0]));
#else
    recv_type[tag][pattern][p] = MPI::FLOAT.Create_indexed(blens.size(), &(blens[0]), &(disp[0]));
#endif
    recv_type[tag][pattern][p].Commit();
  }

  if(verbose)
    cout<<"Got new type\n";
#else
  cerr<<"Error: Not using MPI so types are not created";
#endif
  return 0;
}

int FLComms::GetGEN2UEN(int tag, const int *ENList, int NLocal, const int *gnn2unn,
                        int *gen2uen){
  if(verbose)
    cout<<"int FLComms::GetGEN2UEN(tag="<<tag<<", int *gen2uen)\n";
  
  assert(NLocal>0);
  assert_no_pending_communication(__FILE__, __LINE__);

  // Clear any pre-existing data
  UnregisterHalo(-tag);
  
  size_t nelements = element_ownership[tag].size();
  if(IsParallel()){
#ifdef HAVE_MPI
    // First number all the elements which are owned by the local
    // process.
    size_t npelements=0;
    for(size_t i=0;i<nelements;i++)
      if(element_ownership[tag][i]==(int)MyRank)
        npelements++;

    // *** Just like RegisterHalo(...)
    NOwnedNodes[-tag] = npelements;

    int offset = 0;
    MPI::COMM_WORLD.Scan(&npelements, &offset, 1, MPI::INT, MPI::SUM);
    offset-=npelements;

    // *** Just like RegisterHalo(...)
    this->recv[-tag].resize(NProcs);

    // ***
    vector< vector<int> > unknown_uens(NProcs);
    int cnt=0;
    for(size_t i=0;i<nelements;i++){
      int eowner=element_ownership[tag][i];

      if(eowner<0){ // ie - it's part of a halo that extends beyond this halo level
        continue;
      }

      if(eowner==(int)MyRank){
        gen2uen[i] = offset+cnt;
      }else{
        gen2uen[i] = -1;
        
        for(int j=0;j<NLocal;j++){
          unknown_uens[eowner].push_back(gnn2unn[ENList[i*NLocal+j]-1]);
        }
        this->recv[-tag][eowner].push_back(i);
      }
    }
    
#ifndef NDEBUG
    for(int i=0;i<NProcs;i++){
      // cerr<<i<<": this->recv[-tag][i].size()="<<this->recv[-tag][i].size()<<endl
      //     <<"unknown_uens[i].size()="<<unknown_uens[i].size()<<endl;
      
      assert(this->recv[-tag][i].size()==(unknown_uens[i].size()/NLocal));
    }
#endif

    //
    // Get the UEN for the elements not owned (unknown_uens).
    //

    // Send known/owned UEN's to subordinates.
    assert_no_pending_communication(__FILE__, __LINE__);
    int comm_tag = GetNextCommTag();
    vector<MPI::Request> requests;
    for(int i=0;i<NProcs;i++)
      requests.push_back(MPI::COMM_WORLD.Isend(&(unknown_uens[i][0]), unknown_uens[i].size(), MPI::INT, i, comm_tag));
    
    // Receive elements which need to be identified.
    vector< vector<int> > unknown_elements(NProcs);
    for(int i=0;i<NProcs;i++){
      MPI::Status rstatus;
      MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus); 
      int rcnt=rstatus.Get_count(MPI::INT);
      int source=rstatus.Get_source();
      unknown_elements[source].resize(rcnt);
      requests.push_back(MPI::COMM_WORLD.Irecv(&(unknown_elements[source][0]), rcnt, MPI::INT, source, comm_tag));
    }
    MPI::Request::Waitall(requests.size(), &(requests[0]));
    requests.clear();
    unknown_uens.clear();

    // Identify (UEN) of received elements by first storing the
    // unknowns in a quick lut so they can be quickly searched. Note
    // that an element may be unknown to multiple processes.
    multimap< vector<int>, pair<int, int> > lut;
    vector< vector<int> > returned_uens(NProcs);
    for(int i=0;i<NProcs;i++){
      int nelm = unknown_elements[i].size()/NLocal;
      returned_uens[i].resize(nelm);
      for(int j=0;j<nelm;j++){
        vector<int> elm;
        for(int k=0;k<NLocal;k++)
          elm.push_back(unknown_elements[i][j*NLocal+k]);
        
        lut.insert(pair<vector<int>, pair<int, int> >(elm, pair<int, int>(i, j)));
      }
    }

    // *** Just like RegisterHalo(...)
    this->send[-tag].resize(NProcs);
    for(int i=0;i<NProcs;i++){
      this->send[-tag][i].resize(unknown_elements[i].size()/NLocal);
    }

    // ***
    for(size_t i=0;i<nelements;i++)
      if(element_ownership[tag][i]==(int)MyRank){
        vector<int> elm;
        for(int j=0;j<NLocal;j++)
          elm.push_back(gnn2unn[ENList[i*NLocal+j]-1]);
        
        pair<multimap< vector<int>, pair<int, int> >::iterator, multimap< vector<int>, pair<int, int> >::iterator> range = lut.equal_range(elm);
        for(multimap< vector<int>, pair<int, int> >::const_iterator it=range.first;it!=range.second;++it){
          returned_uens[it->second.first][it->second.second] = gen2uen[i];
          this->send[-tag][it->second.first][it->second.second] = i;
        }
      }

    // Send back uen's
    assert_no_pending_communication(__FILE__, __LINE__);
    comm_tag = GetNextCommTag();
    for(int i=0;i<NProcs;i++)
      requests.push_back(MPI::COMM_WORLD.Isend(&(returned_uens[i][0]), returned_uens[i].size(), MPI::INT, i, comm_tag));

    // Receive elements which have been identified.
    vector< vector<int> > identified_elements(NProcs);
#ifndef NDEBUG
    set<int> hosts_probed;
#endif
    for(int i=0;i<NProcs;i++){
      MPI::Status rstatus;
      MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus); 
      int rcnt=rstatus.Get_count(MPI::INT);
      int source=rstatus.Get_source();
#ifndef NDEBUG
      assert(hosts_probed.count(source)==0);
      hosts_probed.insert(source);
#endif
      identified_elements[source].resize(rcnt);
      requests.push_back(MPI::COMM_WORLD.Irecv(&(identified_elements[source][0]), rcnt, MPI::INT, source, comm_tag));
    }

    MPI::Request::Waitall(requests.size(), &(requests[0]));
    requests.clear();

    // Fill in remaining UEN's
    for(int i=0;i<NProcs;i++){
      //cerr<<i<<": this->recv[-tag][i].size()="<<this->recv[-tag][i].size()<<endl
      //    <<"identified_elements[i].size()="<<identified_elements[i].size()<<endl;
      
      assert(this->recv[-tag][i].size()==identified_elements[i].size());
      for(int j=0;j<(int)this->recv[-tag][i].size();j++){
        assert(gen2uen[this->recv[-tag][i][j]]==-1);
        gen2uen[this->recv[-tag][i][j]] = identified_elements[i][j];
      }
    }

    for(size_t i=0;i<this->send[-tag].size();i++){
      for(size_t j=0;j<this->send[-tag][i].size();j++){
        this->send[-tag][i][j]++;
      }
    }
    for(size_t i=0;i<this->recv[-tag].size();i++){
      for(size_t j=0;j<this->recv[-tag][i].size();j++){
        this->recv[-tag][i][j]++;
      }
    }
#endif
  }else{
    this->send[-tag] = std::deque<std::vector<int> >();
    this->send[-tag].push_back(std::vector<int>());
    this->recv[-tag] = std::deque<std::vector<int> >();
    this->recv[-tag].push_back(std::vector<int>());
    NOwnedNodes[-tag] = nelements;
  
    for(size_t i=0;i<nelements;i++)
      gen2uen[i] = i;
  }

  return 0;
}

int FLComms::GetGNN2UNN(int tag, int NNodes, int NPrivateNodes, int nblocks, int *gnn2unn){
  if(verbose)
    cout<<"int FLComms::GetGNN2UNN(tag="<<tag<<", NNodes="<<NNodes
        <<", NPrivateNodes="<<NPrivateNodes<<", int *gnn2unn)\n";
  
  // This is needed here in the event that this is a serial run.
  Init();
  
  assert_no_pending_communication(__FILE__, __LINE__);
  
#ifdef DEBUG_FLCOMMS
  // Only Initialise array if debugging.
  for(int i=0;i<NNodes*nblocks;i++)
    gnn2unn[i] = -1;
#endif
  
  int d_of_f = NPrivateNodes*nblocks;

  // Get offset value for each partition
  int offset=0;
#ifdef HAVE_MPI
  if(IsParallel()){
    MPI::COMM_WORLD.Scan(&d_of_f, &offset, 1, MPI::INT, MPI::SUM);
    offset-=d_of_f;
  }
#endif
  
  // Number owned nodes
  for(int i=0;i<nblocks;i++){
    for(int j=0; j<NPrivateNodes; j++){
      gnn2unn[i*NNodes+j] = offset + i*NPrivateNodes + j;
    }
  }
  
#ifdef HAVE_MPI
  if(IsParallel()){
    assert(HaveTag(tag));
    
    // Send unn numbering of halo nodes to processes that have these as
    // ghost nodes.
    int comm_tag = GetNextCommTag();
    vector< vector<int> > sendmsg(NProcs);
    vector<MPI::Request> srequests;
    for(size_t i=0; i<(size_t)NProcs; i++){
      for(int j=0; j<nblocks; j++){
        for(vector<int>::const_iterator it=send[tag][i].begin(); it!=send[tag][i].end(); it++){
          // Watch out for Fortran indexing
          int gnn = j*NNodes + (*it) - 1;
          assert(gnn2unn[gnn]!=-1);
          sendmsg[i].push_back(gnn2unn[gnn]);
        }
      }
      srequests.push_back(MPI::COMM_WORLD.Isend(&(sendmsg[i][0]), sendmsg[i].size(), MPI::INT, i, comm_tag));
    }
    
    // Receive and unpack information.
    for(size_t i=0;i<(size_t)NProcs;i++){
      MPI::Status rstatus;
      MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus); 
      int rcnt=rstatus.Get_count(MPI::INT);
      int source=rstatus.Get_source();
      vector<int> recvmsg(rcnt);
      MPI::COMM_WORLD.Recv(&(recvmsg[0]), rcnt, MPI::INT, source, comm_tag);
      
      if(verbose)
        cout<<"assert what? - "<<source<<" :: "<<rcnt<<", "<<nblocks<<", "
            <<recv[tag][source].size()<<", "<<rcnt/nblocks<<" -- "<<recv[1][source].size()
            <<", "<<recv[2][source].size()<<" -- "<<send[1][source].size()<<", "<<send[2][source].size()<<endl;
      
      rcnt/=nblocks;
      assert(recv[tag][source].size()==(size_t)rcnt);
      for(int k=0;k<nblocks;k++){
        for(int j=0;j<rcnt;j++){
          // Watch out for fortran indexing.
          int gnn=k*NNodes+recv[tag][source][j]-1;
          gnn2unn[gnn] = recvmsg[k*rcnt+j];
        }
      }
    }
    
    // Wait for the sends to finish.
    MPI::Request::Waitall(NProcs, &(srequests[0]));
    assert_no_pending_communication(__FILE__, __LINE__);
  }
#endif

  return 0;
}

int FLComms::GetElementOwner(int tag, int eid){
  if(element_ownership.find(tag)==element_ownership.end())
        return MyRank;
  assert(eid<(int)element_ownership[tag].size());
  return element_ownership[tag][eid];
}

int FLComms::GetNextCommTag(){
  fl_comm_tag++;
  if(fl_comm_tag>32767)
    fl_comm_tag = 1;
  return fl_comm_tag;
}

int FLComms::GetNodeOwner(int tag, int nid){
  if(node_ownership.find(tag)==node_ownership.end())
        return MyRank;
  assert(nid<(int)node_ownership[tag].size());
  return node_ownership[tag][nid];
}

int FLComms::GetNOwnedNodes(int tag){
  assert(NOwnedNodes.find(tag)!=NOwnedNodes.end());
  return NOwnedNodes[tag];
}

int FLComms::RegisterHalo(const int& tag, const unsigned int& npnodes, const vector<vector<int> >& send, const vector<vector<int> >& recv){
  if(verbose)
    cout<<"int FLComms::RegisterHalo(const int& tag = "<<tag<<", const unsigned int& npnodes = "<<npnodes<<", const vector<vector<int> >& send, const vector<vector<int> >& recv)\n";
    
  Init();
  
  assert_no_pending_communication(__FILE__, __LINE__);
  
  // Input check
  if(HaveTag(tag)){
    return 0;
  }
  assert((int)send.size()<=NProcs);
  assert((int)recv.size()<=NProcs);
  
  NOwnedNodes[tag]=npnodes;
  
  this->send[tag].resize(NProcs, vector<int>(0));
  this->recv[tag].resize(NProcs, vector<int>(0));
  for(int i=0;i<(int)send.size();i++){
    this->send[tag][i].insert(this->send[tag][i].end(), send[i].begin(), send[i].end());
  }
  for(int i=0;i<(int)recv.size();i++){
    this->recv[tag][i].insert(this->recv[tag][i].end(), recv[i].begin(), recv[i].end());
  }
  
  return 0;
}

int FLComms::RegisterHalo(int tag, int _NOwnedNodes,
                          const int *_send, const int *_bptr_send, 
                          const int *_recv, const int *_bptr_recv){
  if(verbose)
    cout<<"int FLComms::RegisterHalo(tag="<<tag<<" ... )\n";
  
  Init();
  
  assert_no_pending_communication(__FILE__, __LINE__);
  if(HaveTag(tag))
    return 0;
  
  // Work out some basic facts and figures.
  NOwnedNodes[tag] = _NOwnedNodes;

  send[tag].resize(NProcs);
  recv[tag].resize(NProcs);
  
  if(verbose){
    cout<<"bptr_send = ";
    for(size_t i=0;i<=(size_t)NProcs;i++)
      cout<<_bptr_send[i]<<" ";
    cout<<endl;
  }

  for(size_t i=0;i<(size_t)NProcs;i++){
    int slen = _bptr_send[i+1] - _bptr_send[i];
    if(verbose)
      cout<<"send for "<<i<<" has length "<<slen<<endl;

    assert(slen>=0);

    if(slen>0){
      send[tag][i] = vector<int>(slen);
      memcpy(&(send[tag][i][0]), _send+_bptr_send[i]-1, slen*sizeof(int));
    }
    
    int rlen = _bptr_recv[i+1] - _bptr_recv[i];
    if(verbose)
      cout<<"recv for "<<i<<" has length "<<rlen<<endl;

    assert(rlen>=0);

    if(rlen>0){
      recv[tag][i] = vector<int>(rlen);
      memcpy(&(recv[tag][i][0]), _recv+_bptr_recv[i]-1, rlen*sizeof(int));
    }
  }
  
  return 0;
}

int FLComms::UnregisterHalo(const int& tag){
  if(verbose)
    cout<<"int FLComms:UnregiserHalo(tag="<<tag<<")\n";
  
  return Reset(tag);
}

int FLComms::RegisterElements(int tag, int NNodes, int NElements, int NLocal, int *ENList){
  if(verbose)
      cout<<"int FLComms::RegisterElements(tag="<<tag<<" ... )\n";
  
  // A halo for elements is also required. First we must establish
  // ownership. We define the owner of an element to be the process
  // with the minimum rank which owns a node on that element,
  // Create the Unique Node Numbering for linear elements across partitions.
  vector<int> gnn2unn(NNodes, -1);
  GetGNN2UNN(tag, NNodes, NOwnedNodes[tag], 1, &(gnn2unn[0]));
  
  // Establish node ownership. -1 indicates that we're ignoring that
  // node as it probably belongs to an extended halo.
  node_ownership[tag].resize(NNodes);
  fill(node_ownership[tag].begin(), node_ownership[tag].end(), -1);
  for(int i=0;i<NOwnedNodes[tag];i++)
    node_ownership[tag][i] = MyRank;
  
  for(int i=0;i<NProcs;i++)
    for(vector<int>::const_iterator nid=recv[tag][i].begin();nid!=recv[tag][i].end();nid++)
              node_ownership[tag][(*nid)-1] = i;
  
  // Set the minimum node owner of an element as the owner for that
  // element. -1 indicates that we're ignoring that node as it
  // probably belongs to an extended halo.
  element_ownership[tag].resize(NElements);
  for(int i=0;i<NElements;i++){
    size_t nid0 = ENList[i*NLocal] - 1;
    element_ownership[tag][i] = node_ownership[tag][nid0];
    
    for(int j=1;(j<NLocal)&&(element_ownership[tag][i]!=-1);j++){
      size_t nid = ENList[i*NLocal+j] - 1;
      element_ownership[tag][i] = min(element_ownership[tag][i], node_ownership[tag][nid]);
    }
  }
  
  vector<int> gen2uen(NElements);
  GetGEN2UEN(tag, ENList, NLocal, &(gnn2unn[0]), &(gen2uen[0]));
  
  return 0;
}

int FLComms::MergeSurfacePatches(int tag, int NSElements, int NLocal, const int *SENList, int *ids){
#ifdef HAVE_MPI
  // Initialisation is required in the event that this is a serial
  // run.
  Init();

  assert_no_pending_communication(__FILE__, __LINE__);
  
  if(IsParallel()){
    // Offset the numbering of the surface id's
    int max_id=0;
    for(size_t i=0;i<(size_t)NSElements;i++)
      max_id = max(max_id, ids[i]);
    int offset=0;    
    MPI::COMM_WORLD.Scan(&max_id, &offset, 1, MPI::INT, MPI::SUM);
    offset-=max_id;
    for(size_t i=0;i<(size_t)NSElements;i++)
      ids[i] +=offset;
    
    // Node ownership is already known as the volume elements should
    // have already been registered using FLComms::RegisterElements.
    int NNodes = node_ownership[tag].size();
    for(size_t i=0;i<(size_t)NProcs;i++)
      for(vector<int>::const_iterator it=recv[tag][i].begin();it!=recv[tag][i].end();it++)
        NNodes = max(*it, NNodes);
    
    vector<int> gnn2unn(NNodes, -1);
    GetGNN2UNN(tag, NNodes, NOwnedNodes[tag], 1, &(gnn2unn[0]));
    
    for(int idloop=0;idloop<100;idloop++){  
      // Pack what we don't know so we can be enlightened.
      // srecv are the surface element IDs being sent to each process.
      // unknown_ids are the corresponding surface element universal nodes.
      vector< vector<int> > srecv(NProcs);
      vector< vector<int> > unknown_ids(NProcs);
      for(int i=0;i<NSElements;i++){
        // Set of owners of local nodes in this surface element (not including
        // this process)
        set<int> neighs;
        for(int j=0;j<NLocal;j++){
          size_t nid = SENList[i*NLocal+j] - 1;
          if(node_ownership[tag][nid]!=(int)MyRank){
            neighs.insert(node_ownership[tag][nid]);
          }
        }
        // Pack the element universal nodes for surface elements with at least
        // one non-owned node for all processes owning one of the nodes in
        // unknown_ids, and pack the corresponding surface ID in srecv
        for(set<int>::iterator it=neighs.begin();it!=neighs.end();it++){
          for(int j=0;j<NLocal;j++){
            unknown_ids[*it].push_back(gnn2unn[SENList[i*NLocal+j]-1]);
          }
          srecv[*it].push_back(i);
        }
      }
      
      // Send known/owned UEN's to subordinates.
      // Send the surface element universal nodes of surface elements to be
      // communicated
      assert_no_pending_communication(__FILE__, __LINE__);
      int comm_tag = GetNextCommTag();
      vector<MPI::Request> requests;
      for(int i=0;i<NProcs;i++)
        requests.push_back(MPI::COMM_WORLD.Isend(&(unknown_ids[i][0]), unknown_ids[i].size(), MPI::INT, i, comm_tag));
      
      // Receive elements which need to be identified.
      // Receive the surface element universal ndoes of surface elements to be
      // communicated
      vector< vector<int> > unknown_elements(NProcs);
      for(int i=0;i<NProcs;i++){
        MPI::Status rstatus;
        MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus); 
        int rcnt=rstatus.Get_count(MPI::INT);
        int source=rstatus.Get_source();
        unknown_elements[source].resize(rcnt);
        requests.push_back(MPI::COMM_WORLD.Irecv(&(unknown_elements[source][0]), rcnt, MPI::INT, source, comm_tag));
      }
      MPI::Request::Waitall(requests.size(), &(requests[0]));
      requests.clear();
      unknown_ids.clear();
      
      // Identify received received elements by first storing the
      // unknowns in a lut that can be quickly searched. Note
      // that an element may be unknown to multiple processes.
      // Generate a map lut:
      //   (surface element universal nodes) -> (process, surface element no)
      // where the surface element number is its communication index for that
      // process
      multimap< vector<int>, pair<int, int> > lut;
      vector< vector<int> > returned_ids(NProcs);
      if(NLocal>0){
        for(int i=0;i<NProcs;i++){
          int nelm = unknown_elements[i].size()/NLocal;
          returned_ids[i].resize(nelm);
          for(int j=0;j<nelm;j++){
            vector<int> elm;
            for(int k=0;k<NLocal;k++)
              elm.push_back(unknown_elements[i][j*NLocal+k]);
            
            lut.insert(pair<vector<int>, pair<int, int> >(elm, pair<int, int>(i, j)));
          }
        }
      }
      
      // Pack the surface element IDs for communication to other processes in
      // returned_ids (which was allocated of correct size above)
      for(int i=0;i<NSElements;i++){
        vector<int> elm;
        for(int j=0;j<NLocal;j++)
          elm.push_back(gnn2unn[SENList[i*NLocal+j]-1]);
        
        pair<multimap< vector<int>, pair<int, int> >::iterator, multimap< vector<int>, pair<int, int> >::iterator> range = lut.equal_range(elm);
        for(multimap< vector<int>, pair<int, int> >::const_iterator it=range.first;it!=range.second;++it){
          returned_ids[it->second.first][it->second.second] = ids[i];
        }
      }
      
      // Send back ids's
      // Send the surface element IDs to other processes
      assert_no_pending_communication(__FILE__, __LINE__);
      comm_tag = GetNextCommTag();
      for(int i=0;i<NProcs;i++)
        requests.push_back(MPI::COMM_WORLD.Isend(&(returned_ids[i][0]), returned_ids[i].size(), MPI::INT, i, comm_tag));
      
      // Receive elements which have been identified.
      // Receive the surface element IDs from other processes in identified_ids
      vector< vector<int> > identified_ids(NProcs);
      for(int i=0;i<NProcs;i++){
        MPI::Status rstatus;
        MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus); 
        int rcnt=rstatus.Get_count(MPI::INT);
        int source=rstatus.Get_source();
        identified_ids[source].resize(rcnt);
        requests.push_back(MPI::COMM_WORLD.Irecv(&(identified_ids[source][0]), rcnt, MPI::INT, source, comm_tag));
      }
      
      MPI::Request::Waitall(requests.size(), &(requests[0]));
      requests.clear();
      
      // Determine how the patches will be renumbered.
      // Perform the merge
      map<int, int> renumber_patches;
      for(int i=0;i<NProcs;i++){
        for(size_t j=0;j<srecv[i].size();j++){
          if(ids[srecv[i][j]]<=identified_ids[i][j])
            // We already have a lower ID than the incoming ID
            break;

          if(renumber_patches.find(ids[srecv[i][j]])==renumber_patches.end())
            // We have a new ID to map to
            renumber_patches[ids[srecv[i][j]]] = identified_ids[i][j];
          else
            // We have a new ID to map to, but we already had a new ID to map
            // to. Use the lowest incoming ID.
            renumber_patches[ids[srecv[i][j]]] = min(renumber_patches[ids[srecv[i][j]]], identified_ids[i][j]);
        }
      }
      // Remap the surface IDs
      for(int i=0;i<NSElements;i++){
        if(renumber_patches.find(ids[i])!=renumber_patches.end()){
          ids[i] = renumber_patches[ids[i]];
        }
      }
      
      // Count the number of merged surface IDs across all processes
      int update_count=0;
      for(map<int, int>::const_iterator it=renumber_patches.begin();it!=renumber_patches.end();it++)
        if(it->first!=it->second){
          update_count++;
        }
      int pupdate_count;
      MPI::COMM_WORLD.Allreduce(&update_count, &pupdate_count, 1, MPI::INT, MPI::SUM);
      if(pupdate_count==0)
        // Nothing changed in this merge. We have no further merges to perform.
        break;
        
      // We just merged some IDs. This means we have to check for indirect
      // merges. Let's go around again ...
    }
  }
  
  // What happens if you have two process sharing one surface but separated by
  // more than 100 partitions? The surface ID isn't merged ...
  
  assert_no_pending_communication(__FILE__, __LINE__);
#endif
  return 0;
}

int FLComms::Reset(){
  if(verbose)
    cout<<"int FLComms::Reset()\n";
  send.clear();
  recv.clear();

  element_send.clear();
  element_recv.clear();

  NOwnedNodes.clear();
  
  node_ownership.clear();
  element_ownership.clear();

  ClearCache();

  Init();
  return 0;
}

int FLComms::Reset(int tag){
  if(verbose)
    cout<<"int FLComms::Reset(int tag)\n";
  if(!HaveTag(tag))
    return 0;
  
  // Be careful to free the actual MPI types
  send[tag].clear(); send.erase(tag);
  recv[tag].clear(); recv.erase(tag);
 
#ifdef HAVE_MPI
  for(map<pair<int, int>, map<int, MPI::Datatype> >::iterator it=send_type[tag].begin();it!=send_type[tag].end();it++) // over patterns
    for(map<int, MPI::Datatype>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++)                        // over processors
      it2->second.Free();
  
  send_type[tag].clear();
  send_type.erase(tag);
  
  for(map<pair<int, int>, map<int, MPI::Datatype> >::iterator it=recv_type[tag].begin();it!=recv_type[tag].end();it++)     // over patterns
    for(map<int, MPI::Datatype>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++) // over processors
      it2->second.Free();
  recv_type[tag].clear();
  recv_type.erase(tag);
#endif

  return 0;
}

int FLComms::GetNProcs() const{
  return NProcs;
}

int FLComms::SetNProcs(int _NProcs){
  NProcs = _NProcs;
  return 0;
}

int FLComms::Test(int tag, const void *ref, int blk_size, int NFields,
                  int stride, int NNodes, int NPrivate){
  if(verbose)
    cout<<"int FLComms::Test(int tag, const void *ref, int blk_size, int NFields, int stride, int NNodes, int NPrivate)\n";

  // This is needed in the event this is a serial run.
  Init();

  assert_no_pending_communication(__FILE__, __LINE__);
  
  int ierr=0;
#ifdef HAVE_MPI
  if(IsParallel()){
    assert(HaveTag(tag));
    
    int rlen = max(NNodes, max(NNodes,stride)*blk_size*NFields);
#ifdef DOUBLEP
    vector<double> work(rlen, 0.0);
    double *v = (double *)ref;
#else
    vector<float> work(rlen, 0.0);
    float *v = (float *)ref;
#endif
    if(verbose)
      cout<<"Create working array\n";
    
    for(int i=0;i<NFields;i++)
      for(int j=0;j<NPrivate;j++)
        for(int k=0;k<blk_size;k++)
          work[i*stride+j*blk_size+k] = v[i*stride+j*blk_size+k];
    
    // Test volume of information to be sent against the volume to be received
    if(verbose)
      cout<<"Test volume of information to be sent against the volume to be received\n";
    {
      bool passed = true;
      vector<int> scount(NProcs);
      for(size_t i=0;i<(size_t)NProcs;i++){
        scount[i] = send[tag][i].size();
      }
      vector<int> rcount(NProcs, 0);
      MPI::COMM_WORLD.Alltoall(&(scount[0]), 1, MPI::INT, &(rcount[0]), 1, MPI::INT);
      for(size_t i=0;i<(size_t)NProcs;i++){
        if(rcount[i]!=(int)recv[tag][i].size()){
          cerr<<"tag = "<<tag<<", should be receiving "<<recv[tag][i].size()<<" but "<<i<<" only sending "<<rcount[i]<<endl;
          passed = false;
        }
      }
      assert(passed);
    } 
    if(verbose)
      cout<<"...passed.\n"
          <<"Halo communication\n";
    
    Update(tag, &(work[0]), blk_size, NFields, stride);
    
    double min_ele = *min_element(work.begin(), work.end());
    double max_ele = *max_element(work.begin(), work.end());
    double tol=(max_ele-min_ele)*1.0e-6;
    
    for(int i=0;i<NFields;i++){
      for(size_t p=0;p<(size_t)NProcs;p++){
        for(vector<int>::iterator it=recv[tag][p].begin(); it!=recv[tag][p].end();it++){
          assert((*it)>NPrivate);
          assert((*it)<=NNodes);
          int nid = *it - 1;
          for(int k=0;k<blk_size;k++){
            double dx=fabs(work[i*stride+nid*blk_size+k]-v[i*stride+nid*blk_size+k]);
            if(dx>tol){
              for(int l=0;l<blk_size;l++){
                cerr<<":-( "<<MyRank<<" "<<p<<" "<<NPrivate<<" "
                    <<nid<<" "<<NNodes<<" "<<v[i*stride+nid*blk_size+l]<<" <--> "<<work[i*stride+nid*blk_size+l]<<endl;
              }
              ierr=-1;
            }
          }
        }
      }
    }
    
    assert(ierr==0); 
    assert_no_pending_communication(__FILE__, __LINE__);
  }
#endif
  return ierr;
}

int FLComms::Test(int tag, const void *ref, int blk_size, int NFields,
                  int stride, int NNodes, int NPrivate, int NElements, const int *ENList){
  if(verbose)
    cout<<"int FLComms::Test(int tag="<<tag<<", const void *ref, int blk_size="<<blk_size
        <<", int NFields="<<NFields<<", int stride="<<stride<<", int NNodes="<<NNodes
        <<", int NPrivate="<<NPrivate<<", int NElements="<<NElements<<", const int *ENList)\n";

  // This is needed in the event this is a serial run.
  Init();

  assert_no_pending_communication(__FILE__, __LINE__);
  
  int err=0;
#ifdef HAVE_MPI
  if(NProcs>1){    
    assert(Test(tag, ref, blk_size, NFields, stride, NNodes, NPrivate)==0);
    
    // Create the Unique Node Numbering for linear elements across partitions.
    vector<int> gnn2unn(NNodes, -1);
    GetGNN2UNN(tag, NNodes, NPrivate, 1, &(gnn2unn[0]));

    // Cache the ghost nodes into sets for fast lookup
    vector< set<int> > ghost_nodes(NProcs);
    for(size_t i=0; i<(size_t)NProcs; i++)
      for(vector<int>::const_iterator it=recv[tag][i].begin();it!=recv[tag][i].end();it++)
        ghost_nodes[i].insert(*it);

    // Find all the elements that are partly owned by another process.
    deque< vector<int> > sendmsg(NProcs);
    for(int i=0; i<NElements; i++){
      for(size_t j=0; j<(size_t)NProcs; j++){
        for(size_t k=0; k<4; k++){
          if(ghost_nodes[j].find(ENList[i*4+k])!=ghost_nodes[j].end()){
            int n0 = ENList[i*4  ]-1; assert(n0>=0); assert(n0<NNodes);
            int n1 = ENList[i*4+1]-1; assert(n1>=0); assert(n1<NNodes);
            int n2 = ENList[i*4+2]-1; assert(n2>=0); assert(n2<NNodes);
            int n3 = ENList[i*4+3]-1; assert(n3>=0); assert(n3<NNodes);
            
            sendmsg[j].push_back(gnn2unn[n0]);
            sendmsg[j].push_back(gnn2unn[n1]);
            sendmsg[j].push_back(gnn2unn[n2]);
            sendmsg[j].push_back(gnn2unn[n3]);
            break;
          }
        }
      }
    }
    
    for(size_t j=0; j<(size_t)NProcs; j++){
      assert((sendmsg[j].size()%4)==0);
    }

    assert_no_pending_communication(__FILE__, __LINE__);
    
    // Send these elements to the processors that partly own them.
    vector<MPI::Request> srequests;
    int comm_tag = GetNextCommTag();
    for(size_t i=0;i<(size_t)NProcs;i++){
      if(i==MyRank)
        assert(sendmsg[i].size()==0);
      
      srequests.push_back(MPI::COMM_WORLD.Isend(&(sendmsg[i][0]), sendmsg[i].size(), MPI::INT, i, comm_tag));
    }
    
    deque< vector<int> > recvmsg(NProcs);
    for(size_t i=0;i<(size_t)NProcs;i++){
      MPI::Status rstatus;
      MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus);
      int rcnt=rstatus.Get_count(MPI::INT);
      if(verbose)
        cout<<"receive int count = "<<rcnt<<endl;
      assert((rcnt%4)==0);
      int source=rstatus.Get_source();
      recvmsg[source].resize(rcnt);
      MPI::COMM_WORLD.Recv(&(recvmsg[source][0]), rcnt, MPI::INT, source, comm_tag);
    }
    MPI::Request::Waitall(NProcs, &(srequests[0]));
    
    // Free some space.
    sendmsg.clear();
    ghost_nodes.clear();
    
    // Put elements into a set for quick lookup
    set<SimpleTetra4> tetra_set;
    for(int i=0; i<NElements; i++){
      int n0 = ENList[i*4  ]-1; assert(n0>=0); assert(n0<NNodes);
      int n1 = ENList[i*4+1]-1; assert(n1>=0); assert(n1<NNodes);
      int n2 = ENList[i*4+2]-1; assert(n2>=0); assert(n2<NNodes);
      int n3 = ENList[i*4+3]-1; assert(n3>=0); assert(n3<NNodes);
      tetra_set.insert(SimpleTetra4(gnn2unn[n0], gnn2unn[n1], gnn2unn[n2], gnn2unn[n3]));
    }
    
    // Check elements are ok
    assert(SimpleTetra4(0, 1, 2, 3)==SimpleTetra4(1, 2, 3, 0));
    assert(SimpleTetra4(0, 1, 2, 3)!=SimpleTetra4(1, 2, 3, 4));
    assert(SimpleTetra4(0, 1, 2, 3)<SimpleTetra4(1, 2, 3, 4));
    for(size_t i=0; i<(size_t)NProcs; i++){
      for(size_t j=0; j<recvmsg[i].size(); j+=4){
        int n0 = recvmsg[i][j  ];
        int n1 = recvmsg[i][j+1];
        int n2 = recvmsg[i][j+2];
        int n3 = recvmsg[i][j+3];
        
        SimpleTetra4 element(n0, n1, n2, n3);
        if(tetra_set.find(element)==tetra_set.end()){
          cerr<<__FILE__<<", "<<__LINE__<<": halo element is not found - "<<n0<<" "<<n1<<" "<<n2<<" "<<n3<<endl;
          err = -1;
        }
      }
    }

    assert(err==0);
    assert_no_pending_communication(__FILE__, __LINE__);
  }
#endif
  return err;
}

int FLComms::Tetra4ToTetra10(int tagT4, int NPrivateNodes, const int *tetra4, int tagT10, int *tetra10, int NElements){
  if(verbose)
    cout<<"int Tetra4ToTetra10(int, int, const int *, int, int *, int)\n";

  // This is needed in the event this is a serial run.
  Init();
  
  assert_no_pending_communication(__FILE__, __LINE__);
  
#ifdef HAVE_MPI
  if(IsParallel()){
    assert(HaveTag(tagT4));
    if(HaveTag(tagT10)){
      // Replaces assert(!HaveTag(tagT10));
      cerr<<"Warning: Halo tag " << tagT10 << " already assigned - unregistering" << endl;
      UnregisterHalo(tagT10);
    }
#ifndef NDEBUG
    for(int i=0;i<NElements;i++){
      for(size_t j=0;j<10;j++){
        assert(tetra10[i*10+j]>0);
      }
    }
#endif
    
    // Create the set of edges in the overlapping region.
    vector< set< SimpleEdge > > send_edges(NProcs), recv_edges(NProcs);
    map<SimpleEdge, size_t> mid_node_index;
    {
      // Cache the set that is the union of nodes to be sent/received
      // between each process.
      int NNodesT4 = 0;
      for(size_t i=0;i<(size_t)NProcs;i++)
        for(vector<int>::const_iterator it=recv[tagT4][i].begin();it!=recv[tagT4][i].end();it++)
          NNodesT4 = max(*it, NNodesT4);
      
      // Create the Unique Node Numbering for linear elements across partitions.
      vector<int> gnn2unn(NNodesT4, -1);
      GetGNN2UNN(tagT4, NNodesT4, NPrivateNodes, 1, &(gnn2unn[0]));
      
      for(int i=0; i<NElements; i++){
        bool IsHalo=false;
        for(size_t j=0;j<4;j++){
          if(tetra4[i*4+j]>NPrivateNodes){ 
            IsHalo=true;
            break;
          }
        }
      }
      
      // Establish ownership of ghost nodes.
      map<size_t, size_t> owners;
      for(size_t i=0;i<(size_t)NProcs;i++)
        for(vector<int>::const_iterator it=recv[tagT4][i].begin();it!=recv[tagT4][i].end();it++){
          size_t unn = gnn2unn[*it-1];
          assert(owners.count(unn)==0);
          owners[unn] = i;
        }

      // 1. Create mid-node index based of edges - index is local
      //    quadratic element numbering.
      // 2. Partition the edges based on the minimum node owner.
      for(int i=0; i<NElements; i++)
        for(int l0=0; l0<3; l0++)
          for(int l1=l0+1; l1<4; l1++){
            int n0 = gnn2unn[tetra4[i*4+l0]-1];
            assert(n0>=0);

            int n1 = gnn2unn[tetra4[i*4+l1]-1];
            assert(n1>=0);
            
            SimpleEdge edge(n0, n1);
            if(mid_node_index.find(edge)!=mid_node_index.end()){
              if((int)mid_node_index[edge]!=tetra10[i*10 + fl_ilink2(l0+1,l1+1)-1]){
                cerr<<__FILE__<<", "<<__LINE__<<": ERROR - mid_node_index["<<edge<<"]!=tetra10[i*10 + fl_ilink2(l0+1,l1+1)-1])\n";
                exit(-1);
              }
            }else{
              mid_node_index[edge] = tetra10[i*10 + fl_ilink2(l0+1,l1+1)-1];
            }
            
            size_t owner;
            if(owners.count(n0))
              owner = owners[n0];
            else
              owner = MyRank;
            
            if(owners.count(n1))
              owner = min(owner, owners[n1]);
            else
              owner = min(owner, (size_t)MyRank);
            
            if(owner!=MyRank)
              recv_edges[owner].insert(edge);
          }
      
      // Send the list of edges expected to be received from each
      // process, to each of the respective process so that they know
      // what they should be sending.
      vector< vector<unsigned> > sendmsg(NProcs);
      vector<MPI::Request> srequests;
      int comm_tag = GetNextCommTag();
      for(size_t i=0;i<(size_t)NProcs;i++){
        if(i==MyRank){
          assert(recv_edges[i].size()==0);
        }
        for(set<SimpleEdge>::const_iterator it=recv_edges[i].begin(); it!=recv_edges[i].end(); it++){
          sendmsg[i].push_back(it->first);
          sendmsg[i].push_back(it->second);
        }
        srequests.push_back(MPI::COMM_WORLD.Isend(&(sendmsg[i][0]), sendmsg[i].size(), MPI::UNSIGNED, i, comm_tag));
      }

      for(size_t i=0;i<(size_t)NProcs;i++){
        MPI::Status rstatus;
        MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, comm_tag, rstatus);
        int rcnt=rstatus.Get_count(MPI::UNSIGNED);
        assert((rcnt%2)==0);
        int source=rstatus.Get_source();
        vector<unsigned> recvmsg(rcnt);
        MPI::COMM_WORLD.Recv(&(recvmsg[0]), rcnt, MPI::UNSIGNED, source, comm_tag);
        for(int j=0;j<rcnt;j+=2){
          SimpleEdge edge(recvmsg[j], recvmsg[j+1]);
          send_edges[source].insert(edge);
        }
      }
      MPI::Request::Waitall(NProcs, &(srequests[0]));
    }

    // Renumbering from nodes on tetra4 to tetra10
    map<int, int> renum;
    for(int i=0;i<NElements;i++)
      for(int j=0;j<4;j++)
        renum[tetra4[i*4+j]] = tetra10[i*10+fl_ilink2(j+1,j+1)-1];

    // Initalise new halo using verticies and new numbering
    deque< vector<int> > sendT10(NProcs), recvT10(NProcs);
    for(size_t i=0;i<(size_t)NProcs;i++){
      for(vector<int>::const_iterator it=send[tagT4][i].begin();it!=send[tagT4][i].end();it++)
        sendT10[i].push_back(renum[*it]);
      for(vector<int>::const_iterator it=recv[tagT4][i].begin();it!=recv[tagT4][i].end();it++)
        recvT10[i].push_back(renum[*it]);
    }
    
    // Add mid edge nodes
    for(size_t i=0;i<(size_t)NProcs;i++){
      for(set<SimpleEdge>::const_iterator it=recv_edges[i].begin(); it!=recv_edges[i].end(); it++){
        assert(mid_node_index.find(*it)!=mid_node_index.end());
        recvT10[i].push_back(mid_node_index[*it]);
      }
      for(set<SimpleEdge>::const_iterator it=send_edges[i].begin(); it!=send_edges[i].end(); it++){
        assert(mid_node_index.find(*it)!=mid_node_index.end());
        sendT10[i].push_back(mid_node_index[*it]);
      }
    }
    
    //
    // Finally - renumber
    //
    
    // Find elements that are in the second layer of shared elements.
    set<int> elements_ghost2;
    for(int i=0;i<NElements;i++){
      bool h2 = true;
      for(size_t j=0; j<4; j++){
        h2 = !(tetra4[i*4+j]<=NPrivateNodes);
        if(!h2)
          break;
      }
      if(h2)
        elements_ghost2.insert(i);
    }
    
    // Buffer halo nodes from tetra10
    set<int> qhalo;
    for(size_t i=0;i<(size_t)NProcs;i++){
      for(vector<int>::iterator it=recvT10[i].begin();it!=recvT10[i].end();it++)
        qhalo.insert(*it);
    }
    
    // Form sorted sets
    // Split the halo into two different sets to make life easy when debugging
    int NNodesT10 = 0;
    set<int> qhalo2;    
    for(int i=0; i<NElements; i++){
      if(elements_ghost2.count(i))
        for(int j=0; j<10; j++){
          size_t gnn = tetra10[i*10+j];
          if(qhalo.count(gnn)==0){
            qhalo2.insert(gnn);
          }
        }
      for(int j=0; j<10; j++)
        NNodesT10 = max(NNodesT10, tetra10[i*10+j]);
    }
    
    // Form renumbering map
    map<int,int> qrenumbering;
    int NPrivateNodesT10;
    {
      int nid = 1; 
      for(int i=1;i<=NNodesT10;i++){
        if((qhalo.count(i)+qhalo2.count(i))==0)
          qrenumbering[i] = nid++;
      }
      NPrivateNodesT10 = nid-1;
      for(set<int>::iterator it=qhalo.begin();it!=qhalo.end();it++){
        assert(qrenumbering.count(*it)==0);
        qrenumbering[*it] = nid++;
      }
      for(set<int>::iterator it=qhalo2.begin();it!=qhalo2.end();it++){
        assert(qrenumbering.count(*it)==0);
        qrenumbering[*it] = nid++;
      }
      assert(NNodesT10==(nid-1));
    }
    
    // Renumber everything
    for(int i=0;i<NElements;i++)
      for(int j=0;j<10;j++){
        assert(qrenumbering.count(tetra10[i*10+j]));
        tetra10[i*10+j] = qrenumbering[tetra10[i*10+j]];
      }
    for(size_t i=0;i<(size_t)NProcs;i++){
      for(vector<int>::iterator it=recvT10[i].begin();it!=recvT10[i].end();it++){
        assert(qrenumbering.count(*it));
        *it = qrenumbering[*it];
      }
      for(vector<int>::iterator it=sendT10[i].begin();it!=sendT10[i].end();it++){
        assert(qrenumbering.count(*it));
        *it = qrenumbering[*it];
      }
    }
    
    // Add halo
    send[tagT10] = sendT10;
    recv[tagT10] = recvT10;
   
    NOwnedNodes[tagT10] = NPrivateNodesT10;
   
    // Need to prevent collision with other isend/prove/recv's.
    MPI::COMM_WORLD.Barrier();
    
    return NPrivateNodesT10;
  }

#endif
  return 0;
}

//#define ENABLE_MPI_ERROR_CKECKING
int FLComms::Update(int tag, void *buf, int blk_size, int NFields, int stride){
  if(verbose)
    cout<<"int FLComms::Update(int tag="<<tag<<", void *buf, int blk_size="<<blk_size<<", int NFields="<<NFields<<")\n";
  
  // This is needed in the event this is a serial run.
  Init();
  
  assert_no_pending_communication(__FILE__, __LINE__);
  
  int comm_tag = GetNextCommTag();

#ifdef HAVE_MPI
  if(IsParallel()){
    assert(HaveTag(tag));
    assert(buf!=NULL);
    
    CreateTypes(tag, blk_size, NFields, stride);
    
    pair<int, int> pattern(blk_size, NFields);
    
    if(verbose)
      cout<<"Setup non-blocking receives: tag="<<comm_tag<<endl;

    vector<MPI::Request> irequests;

    for(map<int, MPI::Datatype>::const_iterator it=recv_type[tag][pattern].begin(); it!=recv_type[tag][pattern].end(); it++){
      if(verbose)
        cout<<"Receiving message from "<<it->first<<endl;
      irequests.push_back(MPI::COMM_WORLD.Irecv(buf, 1, it->second, it->first, comm_tag));
    }
 
    if(verbose)
      cout<<"Setup non-blocking sends\n";
    
    for(map<int, MPI::Datatype>::const_iterator it=send_type[tag][pattern].begin(); it!=send_type[tag][pattern].end(); it++){
      if(verbose)
        cout<<"Sending message to "<<it->first<<endl;
      irequests.push_back(MPI::COMM_WORLD.Isend(buf, 1, it->second, it->first, comm_tag));
    }
    
    if(verbose)
      cout<<"Wait for all non-blocking communications to finish\n";
    MPI::Request::Waitall(irequests.size(), &(irequests[0]));
  }
#endif
  if(verbose)
    cout<<"updated\n";

  return 0;
}

void FLComms::VerboseOff(){
  verbose = false;
}

void FLComms::VerboseOn(){
  verbose = true;
}

ostream &operator<<(ostream& out, FLComms& comm){
  
  out<<"Start of halo database.\n";
  for(map<int, deque<vector<int> > >::const_iterator it=comm.send.begin();it!=comm.send.end();it++){
    int tag = it->first;
    assert(comm.recv.count(tag));
    out<<"Start of halo, tag "<<tag<<endl;
    
    for(size_t p=0;p<(size_t)comm.NProcs;p++){
      out<<"Sending to "<<p<<":\t";
      for(vector<int>::const_iterator n=comm.send[tag][p].begin(); n!=comm.send[tag][p].end();n++){
        out<<*n<<" ";
      }
      out<<endl;
    }
    
    for(size_t p=0;p<(size_t)comm.NProcs;p++){
      out<<"Receiving from "<<p<<":\t";
      for(vector<int>::const_iterator n=comm.recv[tag][p].begin(); n!=comm.recv[tag][p].end();n++){
        out<<*n<<" ";
      }
      out<<endl;
    }
    out<<"End of halo, tag "<<tag<<endl;
  }
  out<<"End of halo database.\n";
  return out;
}

SimpleEdge::SimpleEdge(size_t a, size_t b){
  first  = min(a,b);
  second = max(a,b);
}

SimpleEdge::SimpleEdge(const SimpleEdge& in){
  *this = in;
}

SimpleEdge::SimpleEdge(const pair<size_t, size_t>& in){
  *this = in;
}

SimpleEdge &SimpleEdge::operator=(const SimpleEdge& in){
  first = in.first;
  second = in.second;
  return *this;
}

SimpleEdge &SimpleEdge::operator=(const pair<size_t, size_t>& in){
  first = min(in.first, in.second);
  second = max(in.first, in.second);
  return *this;
}

bool SimpleEdge::operator==(const SimpleEdge& in) const{
  return (first==in.first)&&(second==in.second);
}

bool SimpleEdge::operator!=(const SimpleEdge& in) const{
  return (first!=in.first)||(second!=in.second);
}

bool SimpleEdge::operator<(const SimpleEdge& in) const{
if(first<in.first)
    return true;

  if(first==in.first)
    if(second<in.second)
      return true;

  return false;
}

ostream &operator<<(ostream& out, const SimpleEdge& in){
  out<<in.first<<" "<<in.second;
  return out;
}

SimpleTetra4::SimpleTetra4(size_t n0, size_t n1, size_t n2, size_t n3){
  assert(n0!=n1); assert(n0!=n2); assert(n0!=n3);
  assert(n1!=n2); assert(n1!=n3);
  assert(n2!=n3);
  
  nodes[0] = n0;
  nodes[1] = n1;
  nodes[2] = n2;
  nodes[3] = n3;
}

SimpleTetra4::SimpleTetra4(const SimpleTetra4 &in){
  *this = in;
}

ostream &operator<<(ostream& out, const SimpleTetra4& in){
  out<<in.nodes[0]<<" "<<in.nodes[1]<<" "<<in.nodes[2]<<" "<<in.nodes[3];
  return out;
}

SimpleTetra4 &SimpleTetra4::operator=(const SimpleTetra4& in){
  for(size_t i=0; i<4; i++)
    nodes[i] = in.nodes[i];
  return *this;
}

bool SimpleTetra4::operator==(const SimpleTetra4& in) const{
  set<int> set0, set1;
  for(size_t i=0;i<4;i++){
    set0.insert(nodes[i]);
    set1.insert(in.nodes[i]);
  }
  
  return set0==set1;
}

bool SimpleTetra4::operator!=(const SimpleTetra4& in) const{
  set<int> set0, set1;
  for(size_t i=0;i<4;i++){
    set0.insert(nodes[i]);
    set1.insert(in.nodes[i]);
  }
  
  return set0!=set1;
}

bool SimpleTetra4::operator<(const SimpleTetra4& in) const{
  set<int> set0, set1;
  for(size_t i=0;i<4;i++){
    set0.insert(nodes[i]);
    set1.insert(in.nodes[i]);
  }
  
  return set0<set1;
}

// 
// Fortran interface.
//

// Filthy dirty global object
FLComms halo_manager;

extern "C" {
#define flcomms_clear_cache_fc F77_FUNC_(flcomms_clear_cache, FLCOMMS_CLEAR_CACHE)
  void flcomms_clear_cache_fc(){
    halo_manager.ClearCache();
  }
  
#define flcomms_halo_registered_fc F77_FUNC_(flcomms_halo_registered, FLCOMMS_HALO_REGISTERED)
  int flcomms_halo_registered_fc(const int* tag){
    return halo_manager.HaveTag(*tag) ? 1 : 0;
  }
  
#define flcomms_get_nhalo_registered_fc F77_FUNC(flcomms_get_nhalo_registered, FLCOMMS_GET_NHALO_REGISTERED)
  int flcomms_get_nhalo_registered_fc() {
    return halo_manager.GetTags().size();
  }
  
#define flcomms_export_halo_fc F77_FUNC(flcomms_export_halo, FLCOMMS_EXPORT_HALO)
  int flcomms_export_halo_fc(const int* tag, int* sends, int* nsends, int* receives, int* nreceives, const int* sends_size, const int* receives_size, const int* nprocs){
    // Assumes sends and receives, and nsends and nreceives, have been allocated
    // to be of sufficient length

    // Export the halo    
    vector<vector<int > > send, recv;
    int ret=halo_manager.ExportHalo(*tag, send, recv);
    if(ret!=0){
      return ret;
    }
    
    // Copy the halo information out of STL data structures
    size_t index=0;
    for(size_t i=0;i<send.size();i++){
      for(size_t j=0;j<send[i].size();j++){
        sends[index++]=send[i][j];
      }
      nsends[i]=send[i].size();
    }
    index=0;
    for(size_t i=0;i<recv.size();i++){
      for(size_t j=0;j<recv[i].size();j++){
        receives[index++]=recv[i][j];
      }
      nreceives[i]=recv[i].size();
    }
    
    return 0;
  } 
  
#define flcomms_get_nprocs_fc F77_FUNC_(flcomms_get_nprocs, FLCOMMS_GET_NPROCS)
  void flcomms_get_nprocs_fc(int* nprocs){
    *nprocs=halo_manager.GetNProcs();
    
    return;
  }
  
#define flcomms_get_info_fc F77_FUNC_(flcomms_get_info, FLCOMMS_GET_INFO)
  void flcomms_get_info_fc(int *tag, int *NProcs, int *num_to_send, int *num_to_recv){
    halo_manager.GetInfo(*tag, NProcs, num_to_send, num_to_recv);
  }
  
#define flcomms_get_gnn2unn_fc F77_FUNC_(flcomms_get_gnn2unn, FLCOMMS_GET_GNN2UNN)
  void flcomms_get_gnn2unn_fc(const int *tag, const int *nnodes,
                              const int *NOwnedNodes, const int *nblocks, int gnn2unn[]){
    halo_manager.GetGNN2UNN(*tag, *nnodes, *NOwnedNodes, *nblocks, gnn2unn);
  }
 
#define flcomms_get_element_owner_fc F77_FUNC_(flcomms_get_element_owner, FLCOMMS_GET_ELEMENT_OWNER)
  void flcomms_get_element_owner_fc(const int *tag, const int *eid, int *owner){
    *owner = halo_manager.GetElementOwner(*tag, (*eid)-1);
  }
  
#define flcomms_get_node_owner_fc F77_FUNC_(flcomms_get_node_owner, FLCOMMS_GET_NODE_OWNER)
  void flcomms_get_node_owner_fc(const int *tag, const int *nid, int *owner){
    *owner = halo_manager.GetNodeOwner(*tag, (*nid-1));
  }

  
#define flcomms_get_nowned_nodes_fc F77_FUNC_(flcomms_get_nowned_nodes, FLCOMMS_GET_NOWNED_NODES)
  void flcomms_get_nowned_nodes_fc(const int *tag, int *nowned_nodes){
    *nowned_nodes = halo_manager.GetNOwnedNodes(*tag);
  }

#define flcomms_merge_surface_patches_fc F77_FUNC_(flcomms_merge_surface_patches, FLCOMMS_MERGE_SURFACE_PATCHES)
  void flcomms_merge_surface_patches_fc(const int *tag, const int *NSElements,
                                        const int *NLocal, const int *SENList, int *ids){
    halo_manager.MergeSurfacePatches(*tag, *NSElements, *NLocal, SENList, ids);
  }
  
#define flcomms_register_halo_fc F77_FUNC_(flcomms_register_halo, FLCOMMS_REGISTER_HALO)
  int flcomms_register_halo_fc(const int* tag, const int* npnodes, const int* sends, const int* sends_size, const int* receives, const int* receives_size, const int* nsends, const int* nreceives, const int* nprocs){
    
    // Check that sends_size and receives_size are consistent with nsends and nreceives
    int send_size_test = 0;
    for(int i=0;i<*nprocs;i++){
      send_size_test+=nsends[i];
    }
    int recv_size_test = 0;
    for(int i=0;i<*nprocs;i++){
      recv_size_test+=nreceives[i];
    }
    assert(*sends_size == send_size_test);
    assert(*receives_size == recv_size_test);
    
    // Copy the halo information into STL data structures
    vector<vector<int > > send(*nprocs), recv(*nprocs);
    for(int i=0;i<*nprocs;i++){
      send[i].resize(nsends[i]);
      recv[i].resize(nreceives[i]);
    }   
    int index=0;
    for(int i=0;i<*nprocs;i++){
      for(int j=0;j<nsends[i];j++){
        send[i][j]=sends[index++];
      }
    }
    index=0;
    for(int i=0;i<*nprocs;i++){
      for(int j=0;j<nreceives[i];j++){
        recv[i][j]=receives[index++];
      }
    }
    
    // Register the halo
    int ret=halo_manager.RegisterHalo(*tag, *npnodes, send, recv);
    if(ret!=0){
      return ret;
    }
  
    return 0;
  }

#define flcomms_unregister_halo_fc F77_FUNC(flcomms_unregister_halo, FLCOMMS_UNREGISTER_HALO)
  int flcomms_unregister_halo_fc(const int* tag){
    return halo_manager.UnregisterHalo(*tag);
  }
  
#define flcomms_register_elements_fc F77_FUNC_(flcomms_register_elements, FLCOMMS_REGISTER_ELEMENTS)
  int flcomms_register_elements_fc(const int *tag, int *NNodes, 
                                    int *NElements, int *NLocal, int *ENList){
    return halo_manager.RegisterElements(*tag, *NNodes, *NElements, *NLocal, ENList);
  }
  
#define flcomms_register_halo_serial_fc F77_FUNC_(flcomms_register_halo_serial, FLCOMMS_REGISTER_HALO_SERIAL)
  void flcomms_register_halo_serial_fc(const int *tag, const int *_NProcs,
                                       const int send[], const int bptr_send[], 
                                       const int recv[], const int bptr_recv[]){
    halo_manager.SetNProcs(*_NProcs);
    halo_manager.RegisterHalo(*tag, 0,
                              send, bptr_send, recv, bptr_recv);
  }
  
#define flcomms_reset_fc F77_FUNC_(flcomms_reset, FLCOMMS_RESET)
  void flcomms_reset_fc(){
    halo_manager.Reset();
  }
  
#define flcomms_reset_tag_fc F77_FUNC_(flcomms_reset_tag, FLCOMMS_RESET_TAG)
  void flcomms_reset_tag_fc(int *tag){
    halo_manager.Reset(*tag);
  }
  
#define flcomms_test_fc F77_FUNC_(flcomms_test, FLCOMMS_TEST)
  void flcomms_test_fc(const int *tag, const void *ref, const int *blk_size, const int *NFields,
                     const int *stride, const int *NNodes, const int *NPrivate,
                     int *ierr){
#ifdef HAVE_MPI
    *ierr = halo_manager.Test(*tag, ref, *blk_size, *NFields,
                              *stride, *NNodes, *NPrivate);
#else
    *ierr = 0;
#endif
  }

#define flcomms_test2_fc F77_FUNC_(flcomms_test2, FLCOMMS_TEST2)
  void flcomms_test2_fc(const int *tag, const void *ref, const int *blk_size, const int *NFields,
                        const int *stride, const int *NNodes, const int *NPrivate, const int *NElements, const int *ENList,
                        int *ierr){
#ifdef HAVE_MPI
    *ierr = halo_manager.Test(*tag, ref, *blk_size, *NFields,
                              *stride, *NNodes, *NPrivate, *NElements, ENList);
#else
    *ierr = 0;
#endif
  }
  
#define flcomms_test_tetra4_to_tetra10_fc F77_FUNC_(flcomms_test_tetra4_to_tetra10, FLCOMMS_TEST_TETRA4_TO_TETRA10)
  void flcomms_test_tetra4_to_tetra10_fc(const int *tag, const int tetra4[], const int tetra10[], 
                                         const int *NElements, const int *NPrivateT10, 
                                         const int *NNodesT10, const void *X, int *ierr){
    *ierr = 0;
#ifdef HAVE_MPI
#ifdef DOUBLEP
    vector<double> ref(*NNodesT10, 0.0);
#else
    vector<float> ref(*NNodesT10, 0.0);
#endif
    for(int i=0;i<*NElements;i++){
      for(int j=0;j<4;j++){
        for(int k=j;k<4;k++){
#ifdef DOUBLEP
          double x0 = ((double *)X)[tetra4[4*i+j]-1];
          double x1 = ((double *)X)[tetra4[4*i+k]-1];
#else
          float x0 = ((float *)X)[tetra4[4*i+j]-1];
          float x1 = ((float *)X)[tetra4[4*i+k]-1];
#endif
          int nid = tetra10[i*10 + fl_ilink2(j+1,k+1)-1]-1;
          assert(nid<(*NNodesT10));
          ref[nid] = (x0+x1)*0.5;
        }
      }
    }
    
    *ierr = halo_manager.Test(*tag, (void *)(&(ref[0])), 1, 1, 1, *NNodesT10, *NPrivateT10);
#endif
  }

#define flcomms_tetra4_to_tetra10_fc F77_FUNC_(flcomms_tetra4_to_tetra10, FLCOMMS_TETRA4_TO_TETRA10)
  void flcomms_tetra4_to_tetra10_fc(int *tagT4, int *NPrivateNodes, const int *tetra4, int *tagT10, int *tetra10, int *NElements, int *NPrivateNodesT10){
#ifdef HAVE_MPI
    *NPrivateNodesT10 = halo_manager.Tetra4ToTetra10(*tagT4, *NPrivateNodes, tetra4, *tagT10, tetra10, *NElements);
#endif
  } 
  
#define flcomms_update_fc F77_FUNC(flcomms_update, FLCOMMS_UPDATE)
  void flcomms_update_fc(int *tag, void *buf, int *blk_size, int *NFields, int *stride){
#ifdef HAVE_MPI
    halo_manager.Update(*tag, buf, *blk_size, *NFields, *stride);
#endif
  }
}
