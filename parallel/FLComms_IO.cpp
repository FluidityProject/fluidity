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

#include "FLComms_IO.h"

using namespace std;

namespace Fluidity{

  HaloReadError ReadHalo(const string& filename, int& process, int& nprocs, map<int, int>& npnodes, map<int, vector<vector<int> > >& send, map<int, vector<vector<int> > >& recv){ 
    // Read the halo file
    TiXmlDocument doc(filename);
    if(!doc.LoadFile()){
      doc.ErrorDesc();
      return HALO_READ_FILE_NOT_FOUND;
    }
    
    const char* charBuffer;
    ostringstream buffer;
     
    // Extract the XML header
    TiXmlNode* header = doc.FirstChild();
    while(header != NULL and header->Type() != TiXmlNode::DECLARATION){
      header = header->NextSibling();
    }

    // Extract the root node
    TiXmlNode* rootNode = header->NextSiblingElement();
    if(rootNode == NULL){
      return HALO_READ_FILE_INVALID;
    }
    TiXmlElement* rootEle = rootNode->ToElement();
    
    // Extract process
    charBuffer = rootEle->Attribute("process");
    if(charBuffer == NULL){
      return HALO_READ_FILE_INVALID;
    }
    process = atoi(charBuffer);
    if(process < 0){
      return HALO_READ_FILE_INVALID;
    }
    
    // Extract nprocs
    charBuffer = rootEle->Attribute("nprocs");
    if(charBuffer == NULL){
      return HALO_READ_FILE_INVALID;
    }
    nprocs = atoi(charBuffer);
    if(process >= nprocs){
      return HALO_READ_FILE_INVALID;
    }
    
    // Extract halo data for each process for each tag
    npnodes.clear();
    send.clear();
    recv.clear();
    // Find the next halo element
    for(TiXmlNode* haloNode = rootEle->FirstChildElement("halo");haloNode != NULL;haloNode = haloNode->NextSiblingElement("halo")){
      if(haloNode == NULL){
        break;
      }
      TiXmlElement* haloEle = haloNode->ToElement();
      
      // Extract the tag
      charBuffer = haloEle->Attribute("tag");
      if(charBuffer == NULL){
        return HALO_READ_FILE_INVALID;
      }
      int tag = atoi(charBuffer);
      send[tag] = vector<vector<int> >(nprocs);
      recv[tag] = vector<vector<int> >(nprocs);

      // Extract n_private_nodes
      charBuffer = haloEle->Attribute("n_private_nodes");
      if(charBuffer == NULL){
        return HALO_READ_FILE_INVALID;
      }
      npnodes[tag] = atoi(charBuffer);
      if(npnodes[tag] < 0){
        return HALO_READ_FILE_INVALID;
      }
      
      // Find the next halo_data element
      for(TiXmlNode* dataNode = haloEle->FirstChildElement("halo_data");dataNode != NULL;dataNode = dataNode->NextSiblingElement("halo_data")){
        if(dataNode == NULL){
          break;
        }
        TiXmlElement* dataEle = dataNode->ToElement();
      
        // Extract the process
        charBuffer = dataEle->Attribute("process");
        if(charBuffer == NULL){
          return HALO_READ_FILE_INVALID;
        }
        int proc = atoi(charBuffer);
        if(proc < 0 or proc >= nprocs){
          return HALO_READ_FILE_INVALID;
        }
        
        // Check that data for this tag and process has not already been extracted
        if(send[tag][proc].size() > 0 or recv[tag][proc].size() > 0){
          return HALO_READ_FILE_INVALID;
        }
        
        // Permit empty send and receive data elements
        send[tag][proc] = vector<int>();
        recv[tag][proc] = vector<int>();
        
        // Extract the send data
        TiXmlNode* sendDataNode = dataEle->FirstChildElement("send");
        if(sendDataNode != NULL){
          TiXmlNode* sendDataTextNode = sendDataNode->FirstChild();
          while(sendDataTextNode != NULL and sendDataTextNode->Type() != TiXmlNode::TEXT){
            sendDataTextNode = sendDataTextNode->NextSibling();
          }
          if(sendDataTextNode != NULL){
            vector<string> tokens;
            Tokenize(sendDataTextNode->ValueStr(), tokens, " ");
            for(size_t i = 0;i < tokens.size();i++){
              send[tag][proc].push_back(atoi(tokens[i].c_str()));
            }
          }
        }
        
        // Extract the receive data
        TiXmlNode* recvDataNode = dataEle->FirstChildElement("receive");
        if(recvDataNode != NULL){
        TiXmlNode* recvDataTextNode = recvDataNode->FirstChild();
          while(recvDataTextNode != NULL and recvDataTextNode->Type() != TiXmlNode::TEXT){
            recvDataTextNode = recvDataTextNode->NextSibling();
          }
          if(recvDataTextNode != NULL){
            vector<string> tokens;
            Tokenize(recvDataTextNode->ValueStr(), tokens, " ");
            for(size_t i = 0;i < tokens.size();i++){
              recv[tag][proc].push_back(atoi(tokens[i].c_str()));
            }
          }
        }
      }
    }

    return HALO_READ_SUCCESS;
  }

  int ReadHalos(const string& basename){
    int process = 0;
    int nprocs = 1;
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      MPI::COMM_WORLD.Barrier();
      process = MPI::COMM_WORLD.Get_rank();
      nprocs = MPI::COMM_WORLD.Get_size();
    }
#endif
    
    int processRead, nprocsRead;
    map<int, int> npnodes;
    map<int, vector<vector<int> > > send, recv;
    
    ostringstream buffer;
    buffer << basename << "_" << process << ".halo";
    HaloReadError ret = ReadHalo(buffer.str(), processRead, nprocsRead, npnodes, send, recv);
    buffer.str("");
    
    int errorCount = 0;
    if(ret == HALO_READ_FILE_NOT_FOUND){
      // Construct an empty halo from the process zero halo file
      // FIXME: Communicate tags from process zero
      buffer << basename << "_" << 0 << ".halo";
      HaloReadError ret = ReadHalo(buffer.str(), processRead, nprocsRead, npnodes, send, recv);
      if(ret == HALO_READ_SUCCESS and processRead == 0 and nprocsRead <= nprocs){
        for(map<int, int>::iterator iter = npnodes.begin();iter != npnodes.end();iter++){
          iter->second = 0;
        }
        for(map<int, vector<vector<int > > >::iterator iter1 = send.begin();iter1 != send.end();iter1++){
          for(vector<vector<int> >::iterator iter2 = iter1->second.begin();iter2 != iter1->second.end();iter2++){
            iter2->clear();
          }
        }
        for(map<int, vector<vector<int > > >::iterator iter1 = recv.begin();iter1 != recv.end();iter1++){
          for(vector<vector<int> >::iterator iter2 = iter1->second.begin();iter2 != iter1->second.end();iter2++){
            iter2->clear();
          }
        }
      }else{
        errorCount++;
      }
    }else if(ret != HALO_READ_SUCCESS){
      errorCount++;
    }else if(processRead != process or nprocsRead > nprocs){
      errorCount++;
    }
    
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      int globalErrorCount = errorCount;
      MPI::COMM_WORLD.Allreduce(&errorCount, &globalErrorCount, 1, MPI::INT, MPI::SUM);
      errorCount = globalErrorCount;
    }
#endif

    if(errorCount == 0){
      // Update the halo manager
      halo_manager.Reset();
      halo_manager.SetNProcs(nprocs);
      for(map<int, vector<vector<int> > >::const_iterator sendIter = send.begin(), recvIter = recv.begin();sendIter != send.end() and recvIter != recv.end();sendIter++, recvIter++){      
        if(halo_manager.RegisterHalo(sendIter->first, npnodes[sendIter->first], sendIter->second, recvIter->second)){
          errorCount++;
#ifdef HAVE_MPI
          if(MPI::Is_initialized()){
            int globalErrorCount = errorCount;
            MPI::COMM_WORLD.Allreduce(&errorCount, &globalErrorCount, 1, MPI::INT, MPI::SUM);
            errorCount = globalErrorCount;
          }
#endif
          break;
        }
      }
    }
    
    return -errorCount;
  }

  int WriteHalo(const string& filename, const unsigned int& process, const unsigned int& nprocs, const map<int, int>& npnodes, const map<int, vector<vector<int> > >& send, const map<int, vector<vector<int> > >& recv){
    // Input check
    assert(process < nprocs);
    assert(send.size() == recv.size());
    for(map<int, vector<vector<int> > >::const_iterator sendIter = send.begin(), recvIter = recv.begin();sendIter != send.end() and recvIter != recv.end(), recvIter != recv.end();sendIter++, recvIter++){
      assert(recv.count(sendIter->first) != 0);
      assert(npnodes.count(sendIter->first) != 0);
      assert(sendIter->second.size() == recvIter->second.size());
    }
    
    TiXmlDocument doc;
    
    ostringstream buffer;
    
    // XML header
    TiXmlDeclaration* header = new TiXmlDeclaration("1.0", "", "");
    doc.LinkEndChild(header);

    // Add root node
    TiXmlElement* rootEle = new TiXmlElement("halos");
    doc.LinkEndChild(rootEle);
    
    // Add process attribute to root node
    buffer << process;
    rootEle->SetAttribute("process", buffer.str());
    buffer.str("");
    
    // Add nprocs attribute to root node
    buffer << nprocs;
    rootEle->SetAttribute("nprocs", buffer.str());
    buffer.str("");
   
    // Add halo data for each tag
    map<int, int>::const_iterator npnodesIter = npnodes.begin();
    for(map<int, vector<vector<int> > >::const_iterator sendTagIter = send.begin(), recvTagIter = recv.begin();sendTagIter != send.end() and recvTagIter != recv.end() and npnodesIter != npnodes.end();sendTagIter++, recvTagIter++, npnodesIter++){
      // Add halo element to root element
      TiXmlElement* haloEle = new TiXmlElement("halo");
      rootEle->LinkEndChild(haloEle);
      
      // Add tag attribute to halo element
      buffer << sendTagIter->first;
      haloEle->SetAttribute("tag", buffer.str());
      buffer.str("");
      
      // Add n_private_nodes attribute to halo element
      buffer << npnodesIter->second;
      haloEle->SetAttribute("n_private_nodes", buffer.str());
      buffer.str("");
      
      // Add halo data for each process for each tag
      int j = 0;
      for(vector<vector<int> >::const_iterator sendProcIter = sendTagIter->second.begin(), recvProcIter = recvTagIter->second.begin();sendProcIter != sendTagIter->second.end() and recvProcIter != recvTagIter->second.end();sendProcIter++, recvProcIter++, j++){
        if(j == (int)nprocs){
          break;
        }

        // Add halo_data element to halo element
        TiXmlElement* dataEle = new TiXmlElement("halo_data");
        haloEle->LinkEndChild(dataEle);
      
        // Add process attribute to data element
        buffer << j;
        dataEle->SetAttribute("process", buffer.str());
        buffer.str("");

        // Add send data to data element
        TiXmlElement* sendDataEle = new TiXmlElement("send");
        dataEle->LinkEndChild(sendDataEle);

        // Add data to send data element
        for(vector<int>::const_iterator sendDataIter = sendProcIter->begin();sendDataIter != sendProcIter->end();sendDataIter++){
          buffer << *sendDataIter << " ";
        }
        TiXmlText* sendData = new TiXmlText(buffer.str());
        sendDataEle->LinkEndChild(sendData);
        buffer.str("");

        // Add receive data to data element
        TiXmlElement* recvDataEle = new TiXmlElement("receive");
        dataEle->LinkEndChild(recvDataEle);
        
        // Add data to receive data element
        for(vector<int>::const_iterator recvDataIter = recvProcIter->begin();recvDataIter != recvProcIter->end();recvDataIter++){
          buffer << *recvDataIter << " ";
        }
        TiXmlText* recvData = new TiXmlText(buffer.str());
        recvDataEle->LinkEndChild(recvData);
        buffer.str("");
      }
    }
    
    return doc.SaveFile(filename) ? 0 : -1;
  }

  int WriteHalos(const string& basename, bool detectTrailing){
    // Export the halo
    int process = 0;
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      process = MPI::COMM_WORLD.Get_rank();
    }
#endif
    size_t nprocs = halo_manager.GetNProcs();
    set<int> haloTags = halo_manager.GetTags();
    map<int, int> npnodes;
    map<int, vector<vector<int> > > send, recv;
    for(set<int>::const_iterator iter = haloTags.begin();iter != haloTags.end();iter++){
      npnodes[*iter] = halo_manager.GetNOwnedNodes(*iter);
      send[*iter] = vector<vector<int> >();
      recv[*iter] = vector<vector<int> >();
      if(halo_manager.ExportHalo(*iter, send[*iter], recv[*iter]) != 0){
        return -1;
      }else if(send[*iter].size() != nprocs or recv[*iter].size() != nprocs){
        return -1;
      }
    }
    
    int nparts;
    if(detectTrailing){
      // Get the number of active partitions - it will be assumed that inactive
      // partitions are all trailing processes
      nparts = 0;
#ifdef HAVE_MPI
      int lnparts;
      if(halo_manager.IsEmpty()){
        lnparts = -1;
      }else{
        lnparts = MPI::COMM_WORLD.Get_rank();  
      }
      MPI::COMM_WORLD.Allreduce(&lnparts, &nparts, 1, MPI::INT, MPI::MAX);
#endif
      nparts++;
#ifdef DDEBUG
      // Assert that inactive partitions are all trailing processes
      if(process < nparts){
        assert(!halo_manager.IsEmpty());
      }else{
        assert(halo_manager.IsEmpty());
      }
#endif
    }else{
      nparts = nprocs;
    }

    // Write out the halos - if trailing partitions are being detected, write
    // for active partitions (only if there is more than one active partition)
    int ret = 0;
    if(not detectTrailing or (nparts > 1 and process < nparts)){
      ostringstream buffer;
      buffer << basename << "_" << process << ".halo";
      ret = WriteHalo(buffer.str(), process, nparts, npnodes, send, recv);
    }
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      int anyWriteFailed = ret ? -1 : 0;
      MPI::COMM_WORLD.Allreduce(&ret, &anyWriteFailed, 1, MPI::INT, MPI::SUM);
      ret = anyWriteFailed;
    }
#endif
    return ret;
  }
  
}

extern "C"{
  int cReadHalos(const char* basename, const int* basename_len){
    return Fluidity::ReadHalos(string(basename, *basename_len));
  }

  int cWriteHalos(const char* basename, const int* basename_len){
    return Fluidity::WriteHalos(string(basename, *basename_len));
  }
}
